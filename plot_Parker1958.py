# -*- coding: utf-8 -*-

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10

import scipy as sp
import scipy.integrate

import matplotlib as mpl
import matplotlib.pyplot as plt

import const



##
## Calculate
##

T_0 = 2e6  ## [K]
print("T_0 = %.2e [K]" % T_0)
## Sound speed
cs = sqrt(2*const.k_B*T_0/const.m_p)
print("cs = %.2e [cm s-1]" % cs)
## Sonic radius
r_s = const.G*const.M_sol/(2*cs**2)
print("r_s = %.2e [cm] = %.2e R_sol" % (r_s, r_s/const.R_sol))


##
## Velocity

def vel_rhs(x, u):
    F = (x-1)/(u**2-1)  if x-1 > 1e-3  else  0.5
    return 2*u/x**2 * F

r = np.linspace(r_s, const.AU, 1000)
solver = sp.integrate.solve_ivp(vel_rhs, t_span=[r[0]/r_s, r[-1]/r_s], y0=[1], t_eval=r/r_s)
v = solver.y[0] * cs


##
## Density

N_0 = 8.8  # [cm^{-3}]
J = r[-1]**2 * N_0 * v[-1]
N = J / (r**2 * v)



##
## Save
##

header = \
""" Parker model for solar wind

 Columns
 -------

 0: r [cm]

 T_0 = 2e6 K
 1: N [cm^{-3}]
 2: v [cm s^{-1}]
"""

arr = np.stack([r, N, v], axis=-1)
np.savetxt("Parker.csv", arr, fmt='%.2e', header=header, comments='#', delimiter=';')



##
## Plot
##

## rc settings (see http://matplotlib.sourceforge.net/users/customizing.html#customizing-matplotlib)
mpl.rc('font', family='serif')
mpl.rc('font', size='8.0')
mpl.rc('text', usetex=True)
mpl.rc('lines', linewidth=0.5)
mpl.rc('axes', linewidth=0.5)
mpl.rc('legend', frameon=False)
mpl.rc('legend', handlelength=2.5)

figwidth = 8.0 / 2.54           ## convert cm to in
figheight = 12.0 / 2.54          ## convert cm to in
mpl.rc('figure', figsize=[figwidth, figheight])

fig, ax = plt.subplots(nrows=3, ncols=1)


ax_ = ax[0]
ax_.loglog(r/const.AU, v/1e5)
ax_.set_xlabel(r"$r$ [AU]")
ax_.set_ylabel(r"$v_\mathrm{w}$ [km s$^{-1}$]")


ax_ = ax[1]
ax_.loglog(r/const.AU, N)
ax_.set_xlabel(r"$r$ [AU]")
ax_.set_ylabel(r"$N_\mathrm{w}$ [cm$^{-3}$]")


ax_ = ax[2]
ax_.semilogy(r/const.AU, N*v**2)
ax_.set_xlabel(r"$r$ [AU]")
ax_.set_ylabel(r"$N_\mathrm{w} v_\mathrm{w}^2$ [cm$^{-1}$ s$^{-2}$]")


plt.tight_layout()
plt.show()
#plt.savefig('Parker1958.pdf')
