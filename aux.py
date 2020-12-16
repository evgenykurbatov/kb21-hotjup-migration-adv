# -*- coding: utf-8 -*-

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10
import scipy as sp
import scipy.optimize

import const



##
## Parameters of the HD209458 system by del Burgo and Allende Prieto (2016)
## https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.1400D
##

## Stellar mass
M_s = 1.148 * const.M_sol  ## (1.148 +/- 0.22) M_sol
## Stellar radius
R_s = 1.2 * const.R_sol  ## (1.2 +/- 0.05) R_sol
## Stellar luminosity
L_s = 1.77 * const.L_sol  ## (1.77 +/- 0.14) L_sol
## Stellar age
#t_age = 3.5e9 * const.yr  ## +/- 1.4e9 yr
t_age = 4e9 * const.yr  ## (3.5 +/- 1.4)*1e9 yr
## Initial mass of the planet
M_p_ini = 0.74 * const.M_jup  ## (0.74 +/- 0.06) M_jup
## Radius of the planet
R_p = 1.41 * const.R_jup  ## (1.41 +/- 0.06) R_jup

## Keplerian orbital frequency
Omega_K = lambda r : sqrt(const.G*M_s/r**3)


##
## Atmosphere loss model by Garcia Munoz (2007)
## https://ui.adsabs.harvard.edu/abs/2007P&SS...55.1426G
## This power law is applicable for a > 0.015 AU.
##

## Current orbit for HD209458b
a_ref = 0.047 * const.AU  ## +/- 0.002 AU
print("P_orb_ref = %g [d]" % (2*pi/Omega_K(a_ref) / (24*3600)))

def dotM_p(t, a, dotM_p_ref):
    ## This power law is applicable for `a > 0.015 AU`:
    dotM = dotM_p_ref * t_age/t * (a_ref/a)**2
    return dotM, -1/t, -2/a


##
## Gas model
##

## Mean molecular weight (assumed to be constant)
mu = 1
## Temperature of gas
T = 1e4  ## K
print("T = %.2e [K]" % T)
## Sound velocity
cs = sqrt(const.RR_gas/mu*T)
print("cs = %.2e [cm/s]" % cs)
## Adiabatic index
gamma = 5/3
print("gamma = %g" % gamma)

## Inner radius of the torus
r_in = lambda M_p, a : ( 1 + (M_p/(3*M_s))**(1/3) ) * a

## `Gravitational radius' by Liffman (2003)
## https://ui.adsabs.harvard.edu/abs/2003PASA...20..337L
r_g = 0.5*(gamma-1)/gamma * const.G*M_s/cs**2
print("r_g = %.2e [cm] = %g AU" % (r_g, r_g/const.AU))



##
## The script is executed as a main program
##

if __name__ == '__main__':

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(nrows=1, ncols=1)


    ax_ = ax
    r = np.linspace(a_ref, 2*const.AU, 500)
    h = cs/Omega_K(r) / r
    ax_.loglog(r/const.AU, h)
    ax_.set_xlabel(r"$r$ [AU]")
    ax_.set_ylabel(r"$H/r$")


    plt.tight_layout()
    plt.show()
