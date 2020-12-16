# -*- coding: utf-8 -*-
##
## No photoevaporation
##

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10
import scipy as sp
import scipy.integrate

import const
from aux import *



path = 'nope'



def dotSigma_pe(t, r, r_0):
    return np.zeros_like(r)
