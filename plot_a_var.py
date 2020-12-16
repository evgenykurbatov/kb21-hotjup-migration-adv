# -*- coding: utf-8 -*-

import sys
import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10

import scipy as sp
import scipy.interpolate

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


def load(dotM_p_ref, a_ini, suffix=""):
    ss = Snapshot()

    fname = '%s/%.0e_%g_%g_%g_%.2f%s.h5' % ('xray', dotM_p_ref, 0.01, 1.5, 1e7, a_ini, suffix)
    with h5py.File(fname, 'r') as f:
        ss.t = f['t'][:]
        ss.a = f['a'][:]

        ss.t = np.append(ss.t, [ss.t[-1]])
        ss.a = np.append(ss.a, [a_ref])

    return ss


p02     = load(5e11, 0.2)
p02_lo  = load(2e11, 0.2)
p02_hi  = load(1e12, 0.2)
p02_min = load(2e11, 0.2, "_min")
p02_max = load(1e12, 0.2, "_max")

p03     = load(5e11, 0.3)
p03_lo  = load(2e11, 0.3)
p03_hi  = load(1e12, 0.3)
p03_min = load(2e11, 0.3, "_min")
p03_max = load(1e12, 0.3, "_max")



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
figheight = 13.0 / 2.54          ## convert cm to in
mpl.rc('figure', figsize=[figwidth, figheight])

fig, ax = plt.subplots(nrows=2, ncols=1)

## [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']


ax_ = ax[0]

## Vertical red bar and the line
rect = mpl.patches.Rectangle((2.1e9, 0), 2.8e9, 1)
ax_.add_collection(PatchCollection([rect], facecolor='k', edgecolor='None', alpha=0.25))
ax_.axvline(3.5e9, ls='-', c='k', alpha=0.5)
ax_.axhline(a_ref/const.AU, ls='-', c='gray')

color = u'#1f77b4'
ax_.fill(np.concatenate((p02_lo.t, p02_hi.t[::-1]))/const.yr,
         np.concatenate((p02_lo.a, p02_hi.a[::-1]))/const.AU,
         facecolor=color, edgecolor=None, alpha=0.25)
ax_.semilogx(p02.t/const.yr, p02.a/const.AU, c=color)

color = u'#ff7f0e'
ax_.fill(np.concatenate((p03_lo.t, p03_hi.t[::-1]))/const.yr,
         np.concatenate((p03_lo.a, p03_hi.a[::-1]))/const.AU,
         facecolor=color, edgecolor=None, alpha=0.25)
ax_.semilogx(p03.t/const.yr, p03.a/const.AU, c=color)

legend = [ mpl.lines.Line2D([], [], color='w', label=r"$\alpha = 0.01$"), \
           mpl.lines.Line2D([], [], color='w', label=r"$\beta = 1.5$") ]
ax_.legend(handles=legend, loc=(0.62, 0.85), handlelength=0)

ax_.set_title(r"Varying $\dot{M}_\mathrm{ref}$ by factor $2$")
ax_.set_xlabel(r"$t$ [yr]")
ax_.set_xlim(1e7, 5e9)
ax_.set_ylabel(r"$a$ [AU]")
ax_.set_ylim(0, 0.34)


ax_ = ax[1]

## Vertical red bar and the line
rect = mpl.patches.Rectangle((2.1e9, 0), 2.8e9, 1)
ax_.add_collection(PatchCollection([rect], facecolor='k', edgecolor='None', alpha=0.25))
ax_.axvline(3.5e9, ls='-', c='k', alpha=0.5)
ax_.axhline(a_ref/const.AU, ls='-', c='gray')

color = u'#1f77b4'
ax_.fill(np.concatenate((p02_min.t, p02_max.t[::-1]))/const.yr,
         np.concatenate((p02_min.a, p02_max.a[::-1]))/const.AU,
         facecolor=color, edgecolor=None, alpha=0.25)
ax_.semilogx(p02.t/const.yr, p02.a/const.AU, c=color)

color = u'#ff7f0e'
ax_.fill(np.concatenate((p03_min.t, p03_max.t[::-1]))/const.yr,
         np.concatenate((p03_min.a, p03_max.a[::-1]))/const.AU,
         facecolor=color, edgecolor=None, alpha=0.25)
ax_.semilogx(p03.t/const.yr, p03.a/const.AU, c=color)

legend = [ mpl.lines.Line2D([], [], color='w', label=r"$\alpha = 0.01$"), \
           mpl.lines.Line2D([], [], color='w', label=r"$\beta = 1.5$") ]
ax_.legend(handles=legend, loc=(0.62, 0.85), handlelength=0)

ax_.set_title(r"Varying $L_\mathrm{X,ref}$ by factor $2$")
ax_.set_xlabel(r"$t$ [yr]")
ax_.set_xlim(1e7, 5e9)
ax_.set_ylabel(r"$a$ [AU]")
ax_.set_ylim(0, 0.34)


plt.tight_layout()
plt.savefig('a_var.pdf')
