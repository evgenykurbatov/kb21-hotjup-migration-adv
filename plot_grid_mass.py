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
    t   = {}
    M_p = {}
    M_t = {}

    for i, _ in enumerate(bp.beta):
        for j, _ in enumerate(bp.alpha):
            for k, _ in enumerate(bp.a_ini):
                t[i,j,k]     = None
                try:
                    fname = '%s/%.0e_%g_%g_%g_%.2f.h5' % (path, dotM_p_ref, bp.alpha[j], bp.beta[i], bp.t_ini, bp.a_ini[k])
                    with h5py.File(fname, 'r') as f:
                        t[i,j,k]   = f['t'][:]
                        M_p[i,j,k] = f['M_p'][:]
                        M_t[i,j,k] = f['M_t'][:]
                except OSError:
                    print("Reading '%s' is failed" % fname)

    return t, M_p, M_t

dotM_p_ref = 5e11
t, M_p, M_t = load('xray', dotM_p_ref)



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

## [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

#dashes = [ [], [6,2], [6,2,1,2], [1,2] ]
dashes = [ [], [8,3], [8,2,1,2], [4,2] ]

fig, ax = plt.subplots(nrows=2, ncols=2)


for i, beta in enumerate(bp.beta):
    for j, alpha in enumerate(bp.alpha):
        ax_ = ax[i,j]

        for k, _ in enumerate(bp.a_ini):
            if t[i,j,k] is not None:
                cond = t[i,j,k]/const.yr > 1.01e7

                #p = ax_.loglog(t[i,j,k][cond]/const.yr, M_t[i,j,k][cond], dashes=dashes[k])
                p = ax_.loglog(t[i,j,k][cond]/const.yr, M_t[i,j,k][cond])

                dM_p = M_p_ini - M_p[i,j,k]
                #ax_.loglog(t[i,j,k][cond]/const.yr, dM_p[cond], c=p[0].get_color(), dashes=dashes[k])
                ax_.loglog(t[i,j,k][cond]/const.yr, dM_p[cond], '--', c=p[0].get_color())

        legend = [ mpl.lines.Line2D([], [], color='w', label=(r"$\alpha = %g$" % alpha)), \
                   mpl.lines.Line2D([], [], color='w', label=(r"$\beta = %g$" % beta)) ]
        ax_.legend(handles=legend, loc=(0.075, 0.42), handlelength=0)
        #ax_.legend(handles=legend, loc='center left', handlelength=0)

        ax_.set_xlim(1e7, 5e9)
        ax_.set_ylim(1e21, 1e29)

        ax__ = ax_.twinx()
        ax__.set_yscale('log')
        ax__.set_ylim(ax_.get_ylim())
        ax__.set_yticklabels([])
        #ax__.set_yticks(ax_.get_yticks())

        #ax_.tick_params(axis='y', which='both', labelleft='off', labelright='off')

        ax_.annotate(r"$\Delta M_\mathrm{p}$", (2e7, 1e28))
        if beta == 1:
            ax_.annotate(r"$M_\mathrm{t}$", (2e8, 3e24))
        else:
            ax_.annotate(r"$M_\mathrm{t}$", (2e8, 1e24))

for ax_ in ax[-1,:]:
    ax_.set_xlabel(r"$t$ [yr]")

for ax_ in ax[:,0]:
    ax_.set_ylabel(r"[g]")


plt.tight_layout()
plt.savefig('grid_mass.pdf')
