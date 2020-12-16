# -*- coding: utf-8 -*-

import sys
import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10
import scipy as sp
import scipy.integrate
import ctypes

import const
from pringle import Pringle



##
## The model
##

class Model(object):
    """
    """

    class Params(ctypes.Structure):
        """Parameters for user-supplied functions of `Pringle` class"""
        _fields_ = [
            ('t',       ctypes.c_double),
            ('a',       ctypes.c_double),
            ('adot',    ctypes.c_double),
            ('r_0',     ctypes.c_double),
            ('Omega_0', ctypes.c_double),
            ('n_0',     ctypes.c_double),
            ('beta',    ctypes.c_double),
            ('m_p',     ctypes.c_double),
            ('dotM_p',  ctypes.c_double * 3),
            ('Sigma_0', ctypes.c_double),
        ]



    def __init__(self, x_, obj=None, **kwargs):
        """"""

        ## Initialize the base model for accretion disk
        self.pringle = Pringle()
        self.pringle.grad_angular_momentum = self.grad_angular_momentum
        self.pringle.viscosity_tensor      = self.viscosity_tensor
        self.pringle.torque                = self.torque
        self.pringle.spec_source           = self.spec_source

        ## Boundary condition
        self.pringle.F_L = 1

        ## Set-up model
        self.pringle.setup(x_)

        ## Positions of the nodes
        self.x_  = self.pringle.r_
        ## Centers of the cells
        self.x   = self.pringle.r
        ## Widths of the cells
        self.dx  = self.pringle.dr
        self.dx_ = self.pringle.dr_


        ##
        ## Set-up parameters by importing from an object or from a dictionary
        ##

        ## Importing from `obj`
        if obj:
            ## Avoid dunders but allow for methods, see https://stackoverflow.com/a/21945171/1038377
            for key in [ _ for _ in dir(obj) if _[:2]+_[-2:] != "____" ]:
                setattr(self, key, getattr(obj, key))
        ## Importing from the `kwargs` dict
        for key in kwargs:
            setattr(self, key, kwargs[key])



    def grad_angular_momentum(self, x):
        """Gradient of the Keplerian specific angular momentum"""
        return 0.5 / sqrt(x)



    def viscosity_tensor(self, x, p):
        """Viscosity tensor (r-phi component) per unit surface density"""
        return - p.n_0 * 1.5 * x**(p.beta -1.5)



    def torque(self, x, p):
        """Tidal torque per unit surface density"""
        xi = p.a/p.r_0
        ## Orbital torque factor, see 1980ApJ...241..425G and 2006RPPh...69..119P
        C_0 = 2.82
        return C_0/pi * p.m_p**2 * xi / (x - xi)**2 * (x**1.5 - xi**1.5) / (x**0.5 - xi**0.5)**3



    def spec_source(self, x, p):
        """This source is caused by the reference frame shift AND by photoevaporation"""

        ## Reference frame shift
        xi = p.r_0/p.a
        Q_rf = - ( p.dotM_p[1] + p.dotM_p[0]/(6*xi*p.m_p**(2/3) * self.M_s) + (p.dotM_p[2] - 0.5/p.a) * p.dota ) / p.Omega_0

        ## Photoevaporation
        Q_pe = - self.kappa_X * self.dotSigma_pe(p.t, p.r_0*x, p.r_0) / p.Omega_0

        return Q_rf + Q_pe



    def dotq(self, _t, q, sigma):
        """R.h.s. for the ODE"""

        ## Stellar age at the moment
        t = self.t_ini + _t
        ## Planetary mass
        M_p   = q[0]
        ## Orbital radius
        a     = q[1]

        x  = self.x
        x_ = self.x_
        dx = self.dx

        r_0 = self.r_in(M_p, a)
        xi = a/r_0

        ##
        ## Atmospheric mass loss

        dotM_p = self.dotM_p(t, a, self.dotM_p_ref)

        ##
        ## Orbit migration

        p = self.Params()
        p.m_p = M_p/self.M_s
        p.a = a
        p.r_0 = r_0
        dota_int = 2*x/sqrt(xi) * self.torque(x, p) * sigma
        dota = sp.integrate.simps(dota_int, x, even='avg')
        dota *= - a * dotM_p[0]/M_p

        return [-dotM_p[0], dota]



    def solve(self, a_ini, a_fin, M_p_ini, t_ini, _t_fin, _t_ss):
        """
        Parameters
        ----------
        a_ini : float
            Initial orbit of the planet [cm].
        a_fin : float
            Final orbit of the planet [cm].
        M_p_ini : float
            Initial mass of the planet [g].
        t_ini : float
            Initial time [s].
        _t_fin : float
            Final time relatively to the initial time [s].
        _t_ss : array_like of floats
            Time points to snapshot relatively to the initial time [s].
        """

        ## Event marker
        def stop_marker(t, q, _):
            return q[1] - a_fin
        stop_marker.terminal  = True
        stop_marker.direction = -1


        ## Initial time
        self.t_ini = t_ini

        ##
        ## Grids

        ## Field's grid (it's expandable)
        q = np.empty((1, 2))
        ## Initial state
        q[0,0]  = M_p_ini
        q[0,1]  = a_ini

        ## Density grid (it's expandable)
        sigma = np.array([])
        ## Initial state
        sigma_old = np.zeros_like(self.x)

        ## Field's grid (it's expandable)
        M_t       = np.array([0])
        sigma_max = np.array([0])
        H_gap     = np.array([0])

        ## Computational time grid (it's expandable)
        ## NB: This time points are intervals from `t_ini`
        _t = np.empty(1)
        ## Initial time point
        _t[0] = 0
        ## Initial time step
        dt_est = 1e3 / self.Omega_K(a_ini)
        print("dt_est = %.2e [yr]" % (dt_est/const.yr))
        dt_last = dt_est
        ## If no constraints, use logarithmic time step
        dlogt = 0.001

        ## Indexes in the '_t' array corresponding to the snapshot times (it's expandable)
        j_ss = np.array([], dtype=np.int)

        ## Index of the current time point
        j = 0
        ## Index of the current snapshot time point
        jt_ss = 0
        ## Flag to make a snapshot
        to_save = False
        ## Flag to stop
        to_stop = False

        ## Solve the model
        print("Compute mass transfer...")

        verbose_every = 500
        while True:

            ##
            ## Pre-processing

            ## Expand grids to the next time point
            _t        = np.append(_t,        [np.nan])
            q         = np.append(q,         [np.empty_like(q[0])], axis=0)
            M_t       = np.append(M_t,       [np.nan])
            sigma_max = np.append(sigma_max, [np.nan])
            H_gap     = np.append(H_gap,     [np.nan])

            ## Get next time point for the current time step
            if (dt_est > 0) and (dt_est < dt_last * 10**dlogt):
                dt = dt_est
            else:
                dt = dt_last * 10**dlogt
            _t[j+1] = _t[j] + dt
            ## Going thru the snapshot time?
            if _t[j+1] >= _t_ss[jt_ss]:
                _t[j+1] = _t_ss[jt_ss]
                to_save = True

            if not (j % verbose_every):
                print("%d: t_ini + %e [yr] = %e [yr] : dt = %.2e [yr]" % (j, _t[j]/const.yr, (t_ini + _t[j])/const.yr, dt/const.yr))


            ##
            ## Mass of the planet and the orbital radius

            ## Advance `q` to the next time point
            sol = sp.integrate.solve_ivp(self.dotq, (_t[j], _t[j+1]), q[j], t_eval=[_t[j], _t[j+1]],
                                         method='BDF', events=stop_marker, args=(sigma_old,),
                                         atol=1e-8, rtol=1e-5, dense_output=True)

            ## If an error occured
            if sol.status == -1:
                print("\tERROR: sol.status=%d, '%s'" % (sol.status, sol.message))
                break

            ## If a termination event occured
            if sol.status == 1:
                ## Set current time to an event time
                _t[j+1] = sol.t_events[0][0]
                ## When not using 'dense_output=True'
                #q[j+1] = sol.y[:,0]
                ## When using 'dense_output=True'
                q[j+1] = sol.sol(_t[j+1])
                print("Event at  t_ini + %e [yr] = %e [yr]"
                      % (_t[j+1]/const.yr, (t_ini + _t[j+1])/const.yr))
                to_stop = True
            else:
                q[j+1] = sol.y[:,1]
            if not (j % verbose_every):
                print("\ta = %g [AU]" % (q[j+1,1]/const.AU))

            ## Are we done?
            if _t[j+1] >= _t_fin:
                _t[j+1] = _t_fin
                q[j+1] = sol.y[:,1]
                print("Finished!")
                to_stop = True


            ##
            ## Density

            ## Preparing auxiliaries for density calculation
            M_p   = q[j+1,0]
            a     = q[j+1,1]
            r_0 = self.r_in(M_p, a)
            Omega_0 = self.Omega_K(r_0)
            p = self.Params()
            p.t = t_ini + _t[j+1]
            p.a = a
            p.r_0 = r_0
            p.Omega_0 = Omega_0
            p.m_p = M_p/self.M_s
            p.dotM_p[:] = self.dotM_p(t_ini + _t[j+1], a, self.dotM_p_ref)
            p.dota = self.dotq(_t[j+1], q[j+1], sigma_old)[1]
            p.n_0 = self.alpha * (self.cs/Omega_0/r_0)**2
            p.beta = self.beta
            p.Sigma_0 = p.dotM_p[0] / (2*pi*r_0**2*Omega_0)

            ## Preparing Pringle
            (A, B, C) = self.pringle.prepare(p)
            ## Advancing `Sigma` to the next time point
            sigma_new = self.pringle.advance(sigma_old, Omega_0*dt, A, B, C, p)

            ## Correcting against negative density values
            sigma_new = np.where(sigma_new > 0, sigma_new, np.zeros_like(sigma_new))


            ##
            ## Torus mass

            x = self.x
            M_t[j+1] = 2*pi*r_0**2 * p.Sigma_0 * sp.integrate.simps(x * sigma_new, x, even='avg')


            ##
            ## Gap width

            sigma_max[j+1] = sigma_new.max()
            ## Gap width at the semi-height level
            i = np.where(sigma_new >= 0.5*sigma_max[j+1])[0][0] - 1
            if i >= 1:
                ## Interpolate gap width as the semi-max position
                H_gap[j+1] = r_0 * (x[i] - 1 + (x[i+1] - x[i]) / (sigma_new[i+1] - sigma_new[i]) * (0.5*sigma_max[j+1] - sigma_new[i]))


            ##
            ## Post-processing

            ## Snapshot the current state
            if to_save or to_stop:
                print("Snapshot %d at  t_ini + %e [yr] = %e [yr]"
                      % (jt_ss, (_t[j+1]/const.yr), (t_ini + _t[j+1])/const.yr))
                ##
                if sigma.size == 0:
                    sigma = np.expand_dims(sigma_new, 0).copy()
                else:
                    sigma = np.append(sigma, [sigma_new], axis=0)
                ## Switching to the next snapshot point
                j_ss = np.append(j_ss, [j+1])
                jt_ss += 1
                to_save = False

            if to_stop:
                break

            ## Estimate the time step size for next iteration
            dt_est = -1

            if not (j % verbose_every):
                ## Characteristic time of the viscosity
                (A, B, C), lam = self.pringle.prepare(p, eigvals=True)
                print("\t%.2e <= |lam| <= %.2e" % (min(np.abs(lam)), max(np.abs(lam))))
                print("\tdt_vi = %.2e [yr]" % (1/(Omega_0*max(np.abs(lam))) / const.yr))
                ## Characteristic time of specific source (rest frame shift and photoevaporation)
                #spec_source = np.max(self.spec_source(self.x, p))
                #print("\tspec_source = %.2e" % spec_source)
                #print("\tdt_spec_source = %.2e [yr]" % (1/(Omega_0*np.abs(spec_source)) / const.yr))
                ## Photoevaporation
                S_pe = np.max(self.kappa_X * self.dotSigma_pe(p.t, p.r_0*x, p.r_0))
                if S_pe:
                    print("\tdt_pe = %.2e [yr]" % (1/S_pe / const.yr))

            ## Switching to the next time point
            j += 1
            dt_last = dt
            sigma_old = sigma_new.copy()


        self._t        = _t
        self.q         = q
        self.sigma     = sigma
        self.j_ss      = j_ss
        self.M_t       = M_t
        self.sigma_max = sigma_max
        self.H_gap     = H_gap
        print("... done.")
