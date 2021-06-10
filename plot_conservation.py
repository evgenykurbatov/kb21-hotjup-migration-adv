# -*- coding: utf-8 -*-

import sys
import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10

import scipy as sp
import scipy.integrate

import h5py

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from aux import *
import batch_params as bp



##
## Read
##

from dataclasses import dataclass

@dataclass
class Snapshot:
    pass


def load(dotM_p_ref, a_ini):
    p = Snapshot()

    fname = '%s/%.0e_%g_%g_%g_%.2f.h5' % ('xray', dotM_p_ref, 0.01, 1.5, 1e7, a_ini)
    with h5py.File(fname, 'r') as f:
        p.t   = f['t'][:]
        p.x   = f['x'][:]
        p.a   = f['a'][:]
        p.r_0 = f['r_0'][:]
        p.M_p = f['M_p'][:]
        p.M_t = f['M_t'][:]
        p.J_t = f['J_t'][:]
        p.dotM_pe = f['dotM_pe'][:]
        p.K_pe    = f['K_pe'][:]
        p.dotJ_pe = f['dotJ_pe'][:]
        p.dotM_b  = f['dotM_b'][:]
        p.K_b     = f['K_b'][:]
        p.dotJ_b  = f['dotJ_b'][:]

        p.J_p = p.M_p * p.a**2 * Omega_K(p.a)

    return p


p = load(5e11, 0.3)



##
## Plot
##

## rc settings (see http://matplotlib.sourceforge.net/users/customizing.html#customizing-matplotlib)
mpl.rc('font', family='serif')
mpl.rc('font', size='6.0')
mpl.rc('text', usetex=True)
mpl.rc('lines', linewidth=0.75)
mpl.rc('axes', linewidth=0.5)
mpl.rc('legend', frameon=False)
mpl.rc('legend', handlelength=2.5)

figwidth  = 8.0 / 2.54           ## convert cm to in
figheight = 13.0 / 2.54          ## convert cm to in
mpl.rc('figure', figsize=[figwidth, figheight])

fig, ax = plt.subplots(nrows=2, ncols=1)

fig.suptitle(r"$a_\mathrm{ini} = 0.3$~AU, $\alpha = 0.01$, $\beta = 1.5$")

## [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']


ax_ = ax[0]

l = ax_.loglog(p.t/const.yr, p.J_p, label="$J_\mathrm{p}$")
ax_.loglog(p.t/const.yr, p.J_p[0] - p.J_p, c=l[0].get_color(), ls='--')

ax_.loglog(p.t/const.yr, p.J_t, label="$J_\mathrm{t}$")

dotJ_p = np.gradient(p.J_p, p.t)
dotJ_t = np.gradient(p.J_t, p.t)

J_b = sp.integrate.cumtrapz(p.dotJ_b, p.t, initial=0.0)
l = ax_.loglog(p.t/const.yr, -J_b, label="$-J_\mathrm{b}$")
ax_.loglog(p.t/const.yr, J_b, c=l[0].get_color(), ls='-.')

J_pe = sp.integrate.cumtrapz(p.dotJ_pe, p.t, initial=0.0)
ax_.loglog(p.t/const.yr, -J_pe, label="$-J_\mathrm{pe}$")

J_tot = p.J_p + p.J_t - J_b - J_pe
ax_.loglog(p.t/const.yr, J_tot, ':k', label="total")

#ax_.text(7e8, 8e49, r"Total", ha='center', va='center')
ax_.text(4e8, 1e50, r"$J_\mathrm{p} + J_\mathrm{t} - J_\mathrm{pe} - J_\mathrm{b}$", ha='center', va='center')
ax_.text(5.5e7, 1.1e49, r"$J_\mathrm{p}$", ha='center', va='center')
ax_.text(1e8, 2e48, r"$\Delta J_\mathrm{p}$", ha='center', va='center')
ax_.text(7e8, 5e44, r"$J_\mathrm{t}$", ha='center', va='center')
ax_.text(5.5e7, 1.3e46, r"$-J_\mathrm{b}$", ha='center', va='center')
ax_.text(7e8, 2e47, r"$-J_\mathrm{pe}$", ha='center', va='center')

ax_.set_xlabel(r"$t$ [yr]")
ax_.set_xlim(xmin=1e7)
ax_.set_ylabel(r"[g cm$^2$ s$^{-1}$]")
ax_.set_yticks([1e42, 1e43, 1e44, 1e45, 1e46, 1e47, 1e48, 1e49, 1e50])
ax_.set_ylim(1e42, 4e50)
#ax_.legend(loc='best')


ax_ = ax[1]

l = ax_.loglog(p.t/const.yr, p.M_p, label="$M_\mathrm{p}$")
ax_.loglog(p.t/const.yr, M_p_ini - p.M_p, c=l[0].get_color(), ls='--')

ax_.loglog(p.t/const.yr, p.M_t, label="$M_\mathrm{t}$")

M_b = sp.integrate.cumtrapz(p.dotM_b, p.t, initial=0.0)
l = ax_.loglog(p.t/const.yr, -M_b, label="$-M_\mathrm{b}$")

M_pe = sp.integrate.cumtrapz(p.dotM_pe, p.t, initial=0.0)
ax_.loglog(p.t/const.yr, -M_pe, label="$-M_\mathrm{pe}$")

M_tot = p.M_p + p.M_t - M_b - M_pe
ax_.loglog(p.t/const.yr, M_tot, ':k', label="total")

ax_.set_xlabel(r"$t$ [yr]")
ax_.set_xlim(xmin=1e7)
ax_.set_ylabel(r"[g]")
ax_.set_yticks([1e22, 1e23, 1e24, 1e25, 1e26, 1e27, 1e28, 1e29, 1e30, 1e31])
ax_.set_ylim(1e22, 2e31)
#ax_.legend(loc='best')

ax_.text(4e8, 4e30, r"$M_\mathrm{p} + M_\mathrm{t} - M_\mathrm{pe} - M_\mathrm{b}$", ha='center', va='center')
ax_.text(5.5e7, 4.1e29, r"$M_\mathrm{p}$", ha='center', va='center')
ax_.text(7e8, 3e28, r"$\Delta M_\mathrm{p}$", ha='center', va='center')
ax_.text(7e8, 1.4e24, r"$M_\mathrm{t}$", ha='center', va='center')
ax_.text(5.5e7, 6e24, r"$-M_\mathrm{b}$", ha='center', va='center')
ax_.text(5.5e7, 7e26, r"$-M_\mathrm{pe}$", ha='center', va='center')


plt.tight_layout(h_pad=3.0)
#plt.show()
plt.savefig('conservation.pdf')
