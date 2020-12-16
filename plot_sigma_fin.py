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
## Read
##

from dataclasses import dataclass

@dataclass
class Snapshot:
    pass


def load(path, dotM_p_ref, alpha, beta, t_ini, a_ini):
    ss = Snapshot()
    ss.alpha = alpha
    ss.beta  = beta

    fname = '%s/%.0e_%g_%g_%g_%.2f.h5' % (path, dotM_p_ref, alpha, beta, t_ini, a_ini)
    with h5py.File(fname, 'r') as f:
        ss.t     = f['t'][:]
        ss.x     = f['x'][:]
        ss.r_0   = f['r_0'][:]
        ss.Sigma = f['Sigma'][:]

    return ss

dotM_p_ref = 5e11
ss_a  = load('xray', dotM_p_ref, 0.001, 1.5, 1e7, 0.2)
ss_b1 = load('xray', dotM_p_ref,  0.01,   1, 1e7, 0.2)
ss_b2 = load('xray', dotM_p_ref,  0.01,   1, 1e7, 0.3)
ss_c1 = load('xray', dotM_p_ref,  0.01, 1.5, 1e7, 0.2)
ss_c2 = load('xray', dotM_p_ref,  0.01, 1.5, 1e7, 0.3)



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

fig, ax = plt.subplots(nrows=1, ncols=1)


ax_ = ax
## [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

ax_.axvline(a_ref/const.AU, ls=':', lw=0.5, c='k')

ss = ss_a
ax_.loglog(ss.r_0[-1]*ss.x / const.AU, ss.Sigma[-1], ls='-', c='#1f77b4',
           label=(r"$\alpha = %g$, $\beta = %g$" % (ss.alpha, ss.beta)))

ss = ss_b1
ax_.loglog(ss.r_0[-1]*ss.x / const.AU, ss.Sigma[-1], ls='-',  c='#ff7f0e',
           label=(r"$\alpha = %g$, \hphantom{$0$}$\beta = %g$" % (ss.alpha, ss.beta)))
ss = ss_b2
ax_.loglog(ss.r_0[-1]*ss.x / const.AU, ss.Sigma[-1], ls='--', c='#ff7f0e')

ss = ss_c1
ax_.loglog(ss.r_0[-1]*ss.x / const.AU, ss.Sigma[-1], ls='-',  c='#2ca02c',
           label=(r"$\alpha = %g$, \hphantom{$0$}$\beta = %g$" % (ss.alpha, ss.beta)))
ss = ss_c2
ax_.loglog(ss.r_0[-1]*ss.x / const.AU, ss.Sigma[-1], ls='--', c='#2ca02c')

ax_.legend()
ax_.set_xlabel(r"$r$ [AU]")
ax_.set_xlim(3e-2, 3e2)
#ax_.set_xlim(4e-2, 8e-2)
ax_.set_ylabel(r"$\Sigma$ [g\,cm$^{-2}$]")
ax_.set_ylim(3e-8, 3e-1)

## Display stroke with annotation
x = np.array([2e-1, 1e0])
y = 3e-5 * (x[0]/x)**1.5
ax_.loglog(x, y, ls='-', c='#ff7f0e')
ax_.annotate(r"$\propto r^{-3/2}$", xy=(sqrt(x[0]*x[1]), 0.6*y[0]),
             ha='center', va='center', rotation=-33)

## Display stroke with annotation
x = np.array([2e-1, 1e0])
y = 3e-6 * (x[0]/x)**2
ax_.loglog(x, y, ls='-', c='#1f77b4')
ax_.loglog(x, 0.4*y, ls='-', c='#2ca02c')
ax_.annotate(r"$\propto r^{-2}$", xy=(sqrt(x[0]*x[1]), 0.4*y[0]),
             ha='center', va='center', rotation=-41)


plt.tight_layout()
plt.savefig('sigma_fin.pdf')
