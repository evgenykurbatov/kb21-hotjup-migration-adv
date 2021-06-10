# -*- coding: utf-8 -*-

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10

import h5py

import const
from aux import *



##
## Main
##

def main(config, dotM_p_ref, a_ini, t_ini, alpha, beta):
    print("%s.main:" % __name__)
    print("\tdotM_p_ref = %.0e [g/s] = %.2e [M_sol/yr]" % (dotM_p_ref, dotM_p_ref/(const.M_sol/const.yr)))
    print("\ta_ini = %g [AU]" % a_ini)
    print("\tt_ini = %g [yr]" % t_ini)
    print("\talpha = %g" % alpha)
    print("\tbeta = %g" % beta)

    ## Prepare parameters for run
    a_ini *= const.AU
    t_ini *= const.yr

    ## Positions of the nodes
    x_ = np.logspace(0, log10(1e4*const.AU/r_in(M_p_ini, a_ini)), 1001)

    from model import Model
    model = Model(x_,
                  obj=config,
                  alpha=alpha, beta=beta, dotM_p_ref=dotM_p_ref, a_ini=a_ini, t_ini=t_ini)
    ## Final orbit
    a_fin = a_ref
    ## Final age
    #_t_fin = t_age - t_ini
    _t_fin = 5e9*const.yr - t_ini
    ## Times to snapshot
    tmp = const.yr * np.array([1e7, 1e8, 1e9, 3e9])
    _t_ss = np.unique(np.concatenate([tmp, [_t_fin]]))
    print("_t_ss [yr] =", _t_ss/const.yr)

    model.solve(a_ini, a_fin, M_p_ini, t_ini, _t_fin, _t_ss)

    _t        = model._t
    x         = model.x
    q         = model.q
    sigma     = model.sigma
    j_ss      = model.j_ss
    M_t       = model.M_t
    J_t       = model.J_t
    dotM_pe   = model.dotM_pe
    K_pe      = model.K_pe
    dotJ_pe   = model.dotJ_pe
    dotM_b    = model.dotM_b
    K_b       = model.K_b
    dotJ_b    = model.dotJ_b
    sigma_max = model.sigma_max
    H_gap     = model.H_gap



    ##
    ## Finalize
    ##

    M_p = q[:,0]
    a   = q[:,1]

    t = _t + t_ini

    r_0 = r_in(M_p, a)
    Omega_0 = Omega_K(r_0)
    Sigma_0 = dotM_p(t, a, dotM_p_ref)[0] / (2*pi*r_0**2*Omega_0)

    ## Make the density distribution dimensional
    Sigma = np.multiply(Sigma_0[j_ss], sigma.T).T
    ## To please h5py
    Sigma = Sigma.astype(np.float64)
    K_pe  = K_pe.astype(np.float64)
    K_b   = K_b.astype(np.float64)

    Sigma_max = Sigma_0 * sigma_max

    r = np.array([ r_0_*x  for r_0_ in r_0 ])



    ##
    ## Write
    ##

    with h5py.File('%s/%.0e_%g_%g_%g_%.2f.h5' % (config.path, dotM_p_ref, alpha, beta, t_ini/const.yr, a_ini/const.AU), 'w') as f:
        ## Star
        f.create_dataset('M_s',         data=M_s) \
         .attrs['comment'] = "Stellar mass [g]"
        f.create_dataset('t_age',       data=t_age) \
         .attrs['comment'] = "Stellar age [s]"
        ## Mass loss by the planet
        f.create_dataset('dotM_p_ref',  data=dotM_p_ref) \
         .attrs['comment'] = "Reference atmospheric mass loss rate [g/s]"
        f.create_dataset('a_ref',       data=a_ref) \
         .attrs['comment'] = "Reference orbit for mass loss rate [cm]"
        ## Gas thermodynamics
        f.create_dataset('mu',          data=mu) \
         .attrs['comment'] = "Mean molecular weight"
        f.create_dataset('T',           data=T) \
         .attrs['comment'] = "Gas temperature [K]"
        f.create_dataset('cs',          data=cs) \
         .attrs['comment'] = "Sound velocity [cm s-1]"
        f.create_dataset('gamma',       data=gamma) \
         .attrs['comment'] = "Adiabatic index"
        ## Turbulence
        f.create_dataset('alpha',       data=alpha) \
         .attrs['comment'] = "Alpha-parameter"
        f.create_dataset('beta',        data=beta) \
         .attrs['comment'] = "Power index for radial dependence of the viscosity coefficient"
        ## Time-dependent quantities
        f.create_dataset('t',           data=t,          compression="gzip") \
         .attrs['comment'] = "Time grid [s]"
        f.create_dataset('x',           data=x,          compression="gzip") \
         .attrs['comment'] = "Radial coordinate grid, in the units of 'r_0'"
        f.create_dataset('M_p',         data=M_p,        compression="gzip") \
         .attrs['comment'] = "Planet mass [g]"
        f.create_dataset('M_t',         data=M_t,        compression="gzip") \
         .attrs['comment'] = "Torus mass [g]"
        f.create_dataset('dotM_pe',     data=dotM_pe,    compression="gzip") \
         .attrs['comment'] = "Total mass loss by photoevaporation [g s-1]"
        f.create_dataset('dotM_b',      data=dotM_b,     compression="gzip") \
         .attrs['comment'] = "Total mass loss by conservative flux [g s-1]"
        f.create_dataset('a',           data=a,          compression="gzip") \
         .attrs['comment'] = "Planetary orbit radius [cm]"
        f.create_dataset('r_0',         data=r_0,        compression="gzip") \
         .attrs['comment'] = "Internal radius of the disk [cm]"
        f.create_dataset('J_t',         data=J_t,        compression="gzip") \
         .attrs['comment'] = "Torus angular momentum [g cm2 s-1]"
        f.create_dataset('dotJ_pe',     data=dotJ_pe,    compression="gzip") \
         .attrs['comment'] = "Total torque by photoevaporation [g cm2 s-2]"
        f.create_dataset('dotJ_b',      data=dotJ_b,     compression="gzip") \
         .attrs['comment'] = "Total torque by conservative flux [g cm2 s-2]"
        f.create_dataset('Sigma_max' ,  data=Sigma_max,  compression="gzip") \
         .attrs['comment'] = "Maximum value of the surface density [g cm-2]"
        f.create_dataset('H_gap',       data=H_gap,      compression="gzip") \
         .attrs['comment'] = "Gap width at the semi-height level [cm]"
        ## Snapshots
        f.create_dataset('j_ss',        data=j_ss) \
         .attrs['comment'] = "Indexes in the 't' array corresponding to the snapshot times"
        f.create_dataset('Sigma',       data=Sigma,      compression="gzip") \
         .attrs['comment'] = "Surface density at the snapshot times [g cm-2]"
        f.create_dataset('K_pe',        data=K_pe,       compression="gzip") \
         .attrs['comment'] = "Cumulative torque distribution by photoevaporation [g cm2 s-2]"
        f.create_dataset('K_b',         data=K_b,        compression="gzip") \
         .attrs['comment'] = "Cumulative torque distribution by conservative flux [g cm2 s-2]"



    ##
    ## Plot
    ##

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(nrows=2, ncols=2)

    ax_ = ax[0,0]
    ax_.loglog(r[-1]/const.AU, Sigma[-1,:])
    ax_.set_xlabel(r"$r$ [AU]")
    ax_.set_ylabel(r"$\Sigma$ [g$/$cm$^2$]")

    ax_ = ax[1,0]
    ax_.loglog(t/const.yr, (M_p_ini - M_p)/const.M_jup, label=r"$\Delta M_\mathrm{p} / M_\mathrm{Jup}$")
    ax_.loglog(t/const.yr, M_t/const.M_jup, label=r"$M_\mathrm{t} / M_\mathrm{Jup}$")
    ax_.set_xlabel(r"$t$ [yr]")
    ax_.legend()

    ax_ = ax[0,1]
    ax_.axhline(a_fin/const.AU, c='k', ls=':')
    ax_.semilogx(t/const.yr, a/const.AU)
    ax_.set_xlabel(r"$t$ [yr]")
    ax_.set_ylabel(r"$a$ [AU]")

    ax_ = ax[1,1]
    ax_.loglog(a/const.AU, H_gap/a)
    ax_.set_xlabel(r"$a$ [AU]")
    ax_.set_ylabel(r"$H_\mathrm{gap}/a$")

    plt.tight_layout()
    plt.savefig('%s/%.0e_%g_%g_%g_%.2f.pdf' % (config.path, dotM_p_ref, alpha, beta, t_ini/const.yr, a_ini/const.AU))


    return 0
