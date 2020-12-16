# -*- coding: utf-8 -*-
##
## https://ui.adsabs.harvard.edu/abs/2012MNRAS.422.1880O
##

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10
import scipy as sp
import scipy.integrate

import matplotlib as mpl
import matplotlib.pyplot as plt

import const



##
##
##

M_s = const.M_sol
L_X = 0.1 * 4*pi*const.AU**2
L_X = 10**27.08
print("L_X = %.2e [erg/s] = %.2e L_sol" % (L_X, L_X/const.L_sol))


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

## `Gravitational radius' by Liffman, https://ui.adsabs.harvard.edu/abs/2003PASA...20..337L
r_g = 0.5*(gamma-1)/gamma * const.G*M_s/cs**2
print("r_g = %.2e [cm] = %g AU" % (r_g, r_g/const.AU))


##
## The model

a = -0.438226
b = -0.10658387
c = 0.5699464
d = 0.010732277
e = -0.131809597
f = -1.32285709

integrand = lambda y : (a*b*exp(b*y) + c*d*exp(d*y) + e*f*exp(f*y)) * exp(-(y/57)**10)
y = np.linspace(0, 500, 5000)
C_y = sp.integrate.simps(integrand(y), y)
print("C_y =", C_y)

C_pe = 4.8e-9 * const.M_sol/const.yr / (2*pi/0.95 * C_y * const.AU**2)
print("C_pe = %.2e [g cm-2 s-1]" % C_pe)


def dotSigma_pe(r, r_0, L_X):
    m = M_s/const.M_sol
    y = 0.95/m * (r - r_0)/const.AU
    C = C_pe * m**(-1.148) * (L_X/1e30)**1.14
    return np.where(r > r_0, C/(r/const.AU) * (a*b*exp(b*y) + c*d*exp(d*y) + e*f*exp(f*y)) * exp(-(y/57)**10), np.zeros_like(r))


def dotM_pe(L_X):
    m = M_s/const.M_sol
    return 4.8e-9 * const.M_sol/const.yr * m**(-0.148) * (L_X/1e30)**1.14



##
## Plot
##

## rc settings (see http://matplotlib.sourceforge.net/users/customizing.html#customizing-matplotlib)
mpl.rc('font', family='serif')
mpl.rc('font', size='8.0')
mpl.rc('text', usetex=True)
mpl.rc('lines', linewidth=0.75)
mpl.rc('axes', linewidth=0.5)
mpl.rc('legend', frameon=False)
mpl.rc('legend', handlelength=2.5)

#figwidth = 8.0 / 2.54           ## convert cm to in
#figheight = 16.0 / 2.54          ## convert cm to in
#mpl.rc('figure', figsize=[figwidth, figheight])

fig, ax = plt.subplots(nrows=1, ncols=1)


ax_ = ax
r = np.logspace(log10(0.045*const.AU), log10(10*const.AU), 1000)
r_0 = 0.1*const.AU
ax_.semilogx(r/const.AU, dotSigma_pe(r, r_0, L_X), ls=':', c='#1f77b4')
ax_.semilogx(r/const.AU, np.where(r > r_g, dotSigma_pe(r, r_0, L_X), np.zeros_like(r)), lw=1.5, c='#1f77b4')
ax_.set_xlabel(r"$r$ [AU]")
ax_.set_ylabel(r"$\dot{\Sigma}$ [g\,cm$^{-2}$\,s$^{-1}$]")
#ax_.legend()


plt.tight_layout()
plt.show()
#plt.savefig('Owen2012_appB2.pdf')
