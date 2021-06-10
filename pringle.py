# -*- coding: utf-8 -*-
"""
Pringle's model of viscous disk
"""

import numpy as np
from numpy import pi, sqrt
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg

from tdma import tdma



class Pringle(object):
    """Pringle's model of viscous disk

    This is a model of Pringle  ``[1]_`` for geometrically thin non-steady
    accretion disk. In addition to the viscous force, as in the original model,
    the angular momentum of external forces is taken into account here.

    .. [1] Pringle J E. Accretion discs in astrophysics // ARA&A, 19:137 (1981).
       <http://adsabs.harvard.edu/abs/1981ARA%26A..19..137P>

    Attributes
    ----------
    r_ : array_like
        Positions of the nodes, $\{r_i\}$.
    r : array_like, shape(len(r_)-1,)
        Centers of the cells, $r_{i+1/2} = 0.5 (r_{i+1} + r_i)$.
    dr : array_like, shape(len(r))
        Widths of the cells, $\Delta_{i+1/2} = r_{i+1} - r_i$.
    dr_ : array_like, shape(len(r_),)
        Widths of the node-centered cells,
        $\Delta_i = 0.5 (\Delta_{i+1/2} + \Delta_{i-1/2}$.
        Note: outside cells are never used, so ``dr_[0] == None`` and
        ``dr_[-1] == None``.
    F_L, F_R : float
        Fluxes at the innermost (``x_[0]``) and outermost (``x_[-1]``)
        boundaries, respectively. If ``None`` (default) a free streaming
        condition is used.

    Methods
    -------
    grad_angular_momentum
        Calculate gradient for specific angular momentum.
    viscosity_tensor
        Calculate the momentum flux density tensor (its $r\phi$-conponent).
    torque
        Calculate the angular momentum of external forces (i.e. torque).
    spec_source
        Calculate the mass specific source which is proportional to density.
    source
        Calculate the mass source.
    """

    ## Basic and auxiliary grids
    r_, r = None, None
    dr, dr_ = None, None

    ## Fluxes at the innermost and autermost boundaries, respectively
    F_L, F_R = None, None



    def setup(self, r_):
        """

        Parameters
        ----------
        r_ : array_like
            Positions of the nodes.
        """

        ## Position of the nodes
        self.r_ = r_

        ## Basic and auxiliary grids
        self.r = 0.5*(self.r_[1:] + self.r_[:-1])
        self.dr = self.r_[1:] - self.r_[:-1]
        self.dr_ = np.concatenate(([None], 0.5*(self.r_[2:] - self.r_[:-2]), [None]))

        ## Auxiliary quantities for differential operator
        r   = self.r
        r_  = self.r_
        dr  = self.dr
        dr_ = self.dr_
        self._phi = 1 / (r * dr)
        vinv = 1.0 / self.grad_angular_momentum(r_)
        self._diffus_m = np.concatenate(([None], vinv[1:-1] * r[:-1]**2/dr_[1:-1], [None]))
        self._diffus_p = np.concatenate(([None], vinv[1:-1] * r[1:]**2/dr_[1:-1],  [None]))
        self._torque   = np.concatenate(([None], vinv[1:-1] * 0.5*r_[1:-1], [None]))



    def prepare(self, p=None, eigvals=False, diagnostic=False):
        """Calculate matrix of differential operator and its spectrum

        Parameters
        ----------
        p : optional
            Parameter for user-supplied functions (e.g. class instance, dict or ctype struct).
        eigvals : bool, optional
            Either to calculate eigenvalues of differential operator or not (default).
        diagnostic : bool, optional
            Either to keep viscosity and torque or not (default).
        """

        r    = self.r
        r_   = self.r_
        dr_  = self.dr_
        _phi  = self._phi
        _diffus_m = self._diffus_m
        _diffus_p = self._diffus_p
        _torque   = self._torque

        W  = self.viscosity_tensor(r, p)
        _T = self.torque(r_, p)

        ## Auxiliary quantities for differential operator
        _psi_m = np.concatenate(([None], _diffus_m[1:-1]*W[:-1] - _torque[1:-1]*_T[1:-1], [None]))
        _psi_p = np.concatenate(([None], _diffus_p[1:-1]*W[1:]  + _torque[1:-1]*_T[1:-1], [None]))

        ## Sub-diagonal
        A = _phi[1:] *  _psi_m[1:-1]
        ## Super-diagonal
        C = _phi[:-1] * _psi_p[1:-1]
        ## Diagonal
        B = - _phi * np.concatenate(([_psi_m[1]], _psi_p[1:-2] + _psi_m[2:-1], [_psi_p[-2]]))

        ## Matrix and its eigenvalues
        if eigvals:
            n = r.size
            L = sp.sparse.diags([A, B, C], [-1, 0, 1], (n, n), dtype=np.float).toarray()
            ret = (A, B, C), sp.linalg.eigvals(L)
        else:
            ret = (A, B, C)

        if diagnostic:
            self.W  = W.copy()
            self._T = _T.copy()

        return ret



    def advance(self, Sigma, dt, A, B, C, p=None):
        """Advance time step

        Parameters
        ----------
        Sigma : ndarray of floats
            Density on current time step.
        dt : float
            Time step size.
        A, B, C : ndarrays of floats
            Diagonals of the matrix of differential operator.
        p : optional
            Parameter for user-supplied functions (e.g. class instance, dict or ctype struct).
        """

        r    = self.r
        r_   = self.r_
        _phi = self._phi
        F_L, F_R = self.F_L, self.F_R

        ## Extremely important to make a copy here,
        ## otherwise it'll be overwritten by `R`
        Sigma_old = Sigma.copy()

        A_ = dt * A
        B_ = dt * B
        C_ = dt * C
        Q = self.spec_source(r, p)
        S = self.source(r, p)

        ## Corrections for boundary conditions
        if F_L is None:
            B_[0] = 0
            C_[0] = 0
            F_L = 0
        if F_R is None:
            A_[-1] = 0
            B_[-1] = 0
            F_R = 0

        B_ += 1 - dt * Q

        ## Right-hand side
        R = Sigma_old + dt * S
        R[0]  += dt * _phi[0] * r_[0]*F_L
        R[-1] += - dt * _phi[-1] * r_[-1]*F_R

        ## Solve using Tridiagonal matrix algorithm,
        ## see <https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm>
        Sigma_new = tdma(A_, B_, C_, R)

        ## Solve using SciPy's routine
        #n = phi.size
        #M = sp.sparse.diags([A_, B_, C_], [-1, 0, 1], (n, n), format='csc', dtype=np.float)
        #Sigma_new = sp.sparse.linalg.spsolve(M, R)

        return Sigma_new



    def diagnostic(self, Sigma):
        """
        """

        r    = self.r
        r_   = self.r_
        dr_  = self.dr_
        _phi  = self._phi
        _diffus_m = self._diffus_m
        _diffus_p = self._diffus_p
        _torque   = self._torque
        W  = self.W
        _T = self._T

        _psi_m = np.concatenate(([None], _diffus_m[1:-1]*W[:-1] - _torque[1:-1]*_T[1:-1], [None]))
        _psi_p = np.concatenate(([None], _diffus_p[1:-1]*W[1:]  + _torque[1:-1]*_T[1:-1], [None]))

        F_L_, F_R_ = self.F_L, self.F_R
        rF_ = np.concatenate(([None], _psi_p[1:-1]*Sigma[1:] - _psi_m[1:-1]*Sigma[:-1], [None]))
        if F_L_ is None:
            rF_[0] = rF_[1]
        else:
            rF_[0] = r_[0]*F_L_
        if F_R_ is None:
            rF_[-1] = rF_[-2]
        else:
            rF_[-1] = r_[-1]*F_R_

        #_r2W = np.concatenate(([None], _r[1:-1]**2 * 0.5*(W[1:]*Sigma[1:] + W[:-1]*Sigma[:-1]), [None]))
        r2W = r**2*W*Sigma

        return rF_, r2W



    def grad_angular_momentum(self, r):
        """"""
        raise NotImplementedError()



    def viscosity_tensor(self, r, p):
        """"""
        raise NotImplementedError()



    def torque(self, r, p):
        """"""
        return np.zeros_like(r)



    def spec_source(self, r, p):
        """"""
        return 0



    def source(self, r, p):
        """"""
        return 0



##
## The source is executed as a main program
##

if __name__ == "__main__":
    ##
    ## Green's function solution to the Pringle's equation

    ## Viscosity
    nu_0 = 1
    beta = 0

    ## Positions of the nodes
    r_0 = 1
    r_ = r_0 * np.logspace(np.log10(0.1), np.log10(2), 1001)
    ## Times' grid
    t_0 = r_0**2/(12*nu_0)
    #t = t_0 * np.array([0, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256])
    t = t_0 * np.array([0, 0.008, 0.032, 0.128, 0.512, 2.048, 8.192, 32.768, 131.072, 524.288])
    ## Surface density grid
    Sigma = np.empty((t.size, r_.size-1))

    ## The model
    model = Pringle()

    ## User-supplied functions
    Pringle.grad_angular_momentum = lambda self, r : 0.5/sqrt(r_)
    Pringle.viscosity_tensor = lambda self, r, _ : - nu_0 * 1.5*r**(beta-1.5)

    ## Set-up model
    model.setup(r_)

    ## Centers of the cells
    r = model.r
    ## Widths of the cells
    dr = model.dr

    ## Initial state
    Sigma[0] = np.zeros_like(r)
    i_0 = np.argmin((r - r_0)**2)
    Sigma[0,i_0] = 1/(2*pi*r[i_0]*dr[i_0])

    ## Boundary flux
    model.F_L = 4

    ## Diagnostic variables
    rF_ = np.empty((t.size, r_.size))
    rF_[0] = np.concatenate(([r_0*model.F_L], np.zeros_like(r)))
    r2W = np.empty((t.size, r_.size-1))
    r2W[0] = np.zeros_like(r)

    ## Solve
    (A, B, C), lam = model.prepare(eigvals=True, diagnostic=True)
    print("min(lam) =", min(lam))
    print("max(lam) =", max(lam))
    for j in range(1, t.size):
        print("t[%d] = %g = %g t_0" % (j, t[j], t[j]/t_0))
        Sigma[j] = model.advance(Sigma[j-1], t[j] - t[j-1], A, B, C)
        rF_[j], r2W[j] = model.diagnostic(Sigma[j])


    ##
    ## Plot

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(nrows=1, ncols=2)

    ax_ = ax[0]
    for j in range(1, t.size):
        ax_.plot(r/r_0, Sigma[j], label=(r"$t/t_0 = %g$" % (t[j]/t_0)))
    ax_.legend(loc='best')
    ax_.set_xlabel(r"$r/r_0$")
    ax_.set_ylabel(r"$\Sigma$")

    ax_ = ax[1]
    for j in range(1, t.size):
        rFLam = 0.5*(rF_[j][:-1] + rF_[j][1:]) * sqrt(r)
        p = ax_.plot(r/r_0, rFLam, ls='-')
        ax_.plot(r/r_0, -r2W[j], c=p[0].get_color(), ls='--')
    ax_.set_xlabel(r"$r/r_0$")
    ax_.set_ylabel(r"solid: $r F r^2 \Omega$, dashed: $-r^2 W$")

    plt.tight_layout()
    plt.show()
