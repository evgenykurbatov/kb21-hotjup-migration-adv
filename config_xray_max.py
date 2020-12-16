# -*- coding: utf-8 -*-
##
## Photoevaporation model
##

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10
import scipy as sp
import scipy.integrate

import const
from aux import *



path = 'xray'



##
## The X-ray luminosity evolution based on two papers:
## i) Louden T, Wheatley P J, Briggs K. Reconstructing the high-energy irradiation of the evaporating hot Jupiter HD 209458b
## https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.2396L
## log10(L_X) = 27.08 +/- 0.07
## ii) Tu L et al. The extreme ultraviolet and X-ray Sun in Time: High-energy evolutionary tracks of a solar-like star
## https://ui.adsabs.harvard.edu/abs/2015A%26A...577L...3T
## L_X ~ t^{-1.42}

L_X = lambda t : 2.0 * 10**(27.08) * (t/t_age)**(-1.42)



##
## The model of photoevaporation
## Owen J E, Clarke C J, Ercolano B, 2012, MNRAS, 422, 1880
## https://ui.adsabs.harvard.edu/abs/2012MNRAS.422.1880O
##

from dataclasses import dataclass

@dataclass
class Photoevaporation:
    pass

pe = Photoevaporation()


def tmpfn(y):
    a = -0.438226
    b = -0.10658387
    c = 0.5699464
    d = 0.010732277
    e = -0.131809597
    f = -1.32285709
    return (a*b*exp(b*y) + c*d*exp(d*y) + e*f*exp(f*y)) * exp(-(y/57)**10)

pe.f = tmpfn

y = np.linspace(0, 500, 5000)
pe.C_y = sp.integrate.simps(tmpfn(y), y, even='avg')
print("pe.C_y =", pe.C_y)

pe.C = 4.8e-9 * const.M_sol/const.yr / (2*pi/0.95 * pe.C_y * const.AU**2)
print("pe.C = %.2e [g cm-2 s-1]" % pe.C)


def dotSigma_pe(t, r, r_0):
    m = M_s/const.M_sol
    y = 0.95/m * (r - r_0)/const.AU
    cond = (r >= r_0) & (r >= r_g)
    return pe.C * m**(-1.148) * (L_X(t)/1e30)**1.14 / (r/const.AU) \
        * np.where(cond, pe.f(y), np.zeros_like(r))


def dotM_pe(t, r_0):
    m = M_s/const.M_sol
    y_min = np.where(r_g >= r_0, 0.95/m * (r_g - r_0)/const.AU, 0)
    y = np.linspace(y_min, 500, 5000).T
    integral = sp.integrate.simps(pe.f(y), y, even='avg')
    return 4.8e-9 * const.M_sol/const.yr * m**(-0.148) * (L_X(t)/1e30)**1.14 * integral/pe.C_y


##
## Ride S K, Walker A B C Jr. Absorption of X-rays in the interstellar medium
## https://ui.adsabs.harvard.edu/abs/1977A%26A....61..339R
##
## Opacity for the solar sbundance for ~ 1 keV photons:
kappa_X = 2e-22 / const.m_p
print("kappa_X = %.2e [cm^2/g]" % kappa_X)
