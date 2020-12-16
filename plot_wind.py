# -*- coding: utf-8 -*-

import sys
import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10

import h5py

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from aux import *



##
## Read simulation data
##

from dataclasses import dataclass

@dataclass
class Snapshot:
    pass


def load(path, dotM_p_ref, alpha, beta, t_ini, a_ini):
    ss = Snapshot()
    ss.alpha = alpha
    ss.beta  = beta
    ss.a_ini = a_ini

    fname = '%s/%.0e_%g_%g_%g_%.2f.h5' % (path, dotM_p_ref, alpha, beta, t_ini, a_ini)
    with h5py.File(fname, 'r') as f:
        ss.M_s       = f['M_s'][()]
        ss.mu        = f['mu'][()]
        ss.cs        = f['cs'][()]
        ss.t         = f['t'][:]
        ss.a         = f['a'][:]
        ss.r_0       = f['r_0'][:]
        ss.Sigma_max = f['Sigma_max'][:]

    ss.H = ss.cs / sqrt(const.G*ss.M_s/ss.r_0**3)
    ss.N = ss.Sigma_max / (2*ss.H * ss.mu*const.m_H)

    return ss


ss = [ load('xray', dotM_p_ref=5e11,  alpha=0.01, beta=1.5, t_ini=1e7, a_ini=a_ini)
       for a_ini in [0.2, 0.3, 0.4, 0.5] ]

kappa_X = 2e-22 / const.m_p



##
## Read Parker model
##

parker_data = np.loadtxt('Parker1958.csv', delimiter=';')
c = 12
r   = parker_data[:,c]
N_w = parker_data[:,c+1]
v_w = parker_data[:,c+2]



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

figwidth = 8.0 / 2.54           ## convert cm to in
figheight = 6.0 / 2.54          ## convert cm to in
mpl.rc('figure', figsize=[figwidth, figheight])

## [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

#dashes = [ [], [6,2], [6,2,1,2], [1,2] ]
dashes = [ [], [8,3], [8,2,1,2], [4,2] ]

fig, ax = plt.subplots(nrows=1, ncols=1)

fig.suptitle(r"$\alpha = 0.01$, $\beta = 1.5$")


ax_ = ax
for i, ss_ in enumerate(ss):
    cond = ss_.t/const.yr > 1.001e7
    #ax_.plot(ss_.r_0[cond]/const.AU, (ss_.N*cs**2)[cond], dashes=dashes[i],
    #         label=(r"$a_\mathrm{ini} = %g$ AU" % ss_.a_ini))
    ax_.plot(ss_.r_0[cond]/const.AU, (const.m_H*ss_.N*cs**2 / (1e-3*const.mbar))[cond])
#ax_.plot(r/const.AU, N_w*v_w**2, '-k', label=r"Parker wind")
ax_.plot(r/const.AU, const.m_p*N_w*v_w**2 / (1e-3*const.mbar), '-k')
ax_.legend()
ax_.set_xscale('log')
ax_.set_xlim(0.047, 1.0)
ax_.set_xlabel(r"$a$ [AU]")
ax_.set_yscale('log')
#ax_.set_ylabel(r"[g cm$^{-1}$ s$^{-2}$]")
ax_.set_ylabel(r"[$\mu$bar]")


plt.tight_layout()
plt.savefig('wind.pdf')
