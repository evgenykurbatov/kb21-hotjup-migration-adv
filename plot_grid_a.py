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
    t = {}
    a = {}

    for i, _ in enumerate(bp.beta):
        for j, _ in enumerate(bp.alpha):
            for k, _ in enumerate(bp.a_ini):
                t[i,j,k] = None
                a[i,j,k] = None
                try:
                    fname = '%s/%.0e_%g_%g_%g_%.2f.h5' % (path, dotM_p_ref, bp.alpha[j], bp.beta[i], bp.t_ini, bp.a_ini[k])
                    with h5py.File(fname, 'r') as f:
                        t[i,j,k] = f['t'][:]
                        a[i,j,k] = f['a'][:]
                except OSError:
                    print("Reading '%s' is failed" % fname)

    return t, a


dotM_p_ref = 5e11
nope_t, nope_a = load('nope', dotM_p_ref)
xray_t, xray_a = load('xray', dotM_p_ref)



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

## [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']


for i, beta in enumerate(bp.beta):
    for j, alpha in enumerate(bp.alpha):
        ax_ = ax[j,i]

        ## Vertical red bar and the line
        rect = mpl.patches.Rectangle((2.1e9, 0), 2.8e9, 1)
        #ax_.add_collection(PatchCollection([rect], facecolor='r', edgecolor='None', alpha=0.25))
        #ax_.axvline(3.5e9, ls='-', c='r', alpha=0.5)
        ax_.add_collection(PatchCollection([rect], facecolor='k', edgecolor='None', alpha=0.25))
        ax_.axvline(3.5e9, ls='-', c='k', alpha=0.5)
        ax_.axhline(a_ref/const.AU, ls='-', c='gray')

        for k, _ in enumerate(bp.a_ini):
            """
            if nope_t[i,j,k] is not None:
                ax_.semilogx(nope_t[i,j,k]/const.yr, nope_a[i,j,k]/const.AU, ':', c='#1f77b4')
            if xray_t[i,j,k] is not None:
                ax_.semilogx(xray_t[i,j,k]/const.yr, xray_a[i,j,k]/const.AU, '-', c='#1f77b4')
            """
            if xray_t[i,j,k] is not None:
                p = ax_.semilogx(xray_t[i,j,k]/const.yr, xray_a[i,j,k]/const.AU, '-')
                ax_.semilogx(nope_t[i,j,k]/const.yr, nope_a[i,j,k]/const.AU, ':', c=p[0].get_color())

        legend = [ mpl.lines.Line2D([], [], color='w', label=(r"$\alpha = %g$" % alpha)), \
                   mpl.lines.Line2D([], [], color='w', label=(r"$\beta = %g$" % beta)) ]
        ax_.legend(handles=legend, loc=(0.62, 0.85), handlelength=0)

        ax_.set_xlim(bp.t_ini, 5e9)
        ax_.set_ylim(0, 0.599)

for ax_ in ax[-1,:]:
    ax_.set_xlabel(r"$t$ [yr]")

for ax_ in ax[:,0]:
    ax_.set_ylabel(r"$a$ [AU]")


plt.tight_layout()
plt.savefig('grid_a.pdf')
