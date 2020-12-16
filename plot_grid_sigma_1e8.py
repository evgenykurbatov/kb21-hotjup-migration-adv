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
import batch_params as bp



##
## Read
##

def load(path, dotM_p_ref):
    t     = {}
    x     = {}
    r_0   = {}
    j_ss  = {}
    Sigma = {}

    for i, _ in enumerate(bp.beta):
        for j, _ in enumerate(bp.alpha):
            for k, _ in enumerate(bp.a_ini):
                t[i,j,k]     = None
                try:
                    fname = '%s/%.0e_%g_%g_%g_%.2f.h5' % (path, dotM_p_ref, bp.alpha[j], bp.beta[i], bp.t_ini, bp.a_ini[k])
                    with h5py.File(fname, 'r') as f:
                        t[i,j,k]     = f['t'][:]
                        x[i,j,k]     = f['x'][:]
                        r_0[i,j,k]   = f['r_0'][:]
                        j_ss[i,j,k]  = f['j_ss'][:]
                        Sigma[i,j,k] = f['Sigma'][:]
                except OSError:
                    print("Reading '%s' is failed" % fname)

    return t, x, r_0, j_ss, Sigma

dotM_p_ref = 5e11
t, x, r_0, j_ss, Sigma = load('xray', dotM_p_ref)



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

figwidth = 16.0 / 2.54           ## convert cm to in
figheight = 12.0 / 2.54          ## convert cm to in
mpl.rc('figure', figsize=[figwidth, figheight])

fig, ax = plt.subplots(nrows=2, ncols=2)


for i, beta in enumerate(bp.beta):
    for j, alpha in enumerate(bp.alpha):
        ax_ = ax[i,j]

        for k, _ in enumerate(bp.a_ini):
            ## [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
            if t[i,j,k] is not None:
                n_ss = 1
                n = j_ss[i,j,k][n_ss]
                r = r_0[i,j,k][n]*x[i,j,k]
                #ax_.loglog(r/const.AU, Sigma[i,j,k][n_ss], c='#1f77b4')
                ax_.loglog(r/const.AU, Sigma[i,j,k][n_ss])

        ax_.text(1.6e-1, 3e-4, r"$0.2$ AU", ha='center', va='center')
        ax_.text(9e-1, 4e-7, r"$0.5$ AU", ha='center', va='center')

        legend = [ mpl.lines.Line2D([], [], color='w', label=(r"$\alpha = %g$" % alpha)), \
                   mpl.lines.Line2D([], [], color='w', label=(r"$\beta = %g$" % beta)) ]
        ax_.legend(handles=legend, handlelength=0)

        ax_.set_xlim(1e-1, 3e1)
        ax_.set_ylim(7e-8, 1e-3)

for ax_ in ax[-1,:]:
    ax_.set_xlabel(r"$r$ [AU]")

for ax_ in ax[:,0]:
    ax_.set_ylabel(r"$\Sigma$ [g\,cm$^{-2}$]")


plt.tight_layout()
plt.savefig('grid_sigma_1e8.pdf')
