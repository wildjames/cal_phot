#!/usr/bin/env python
from hipercam.hlog import Hlog

from astropy import coordinates as coord, units as u
from astropy.time import Time
from astropy.convolution import Box1DKernel, convolve
from astropy.stats import sigma_clipped_stats

from scipy.signal import medfilt
from scipy.optimize import minimize, leastsq

import celerite
from celerite import terms
from celerite.modeling import Model

import numpy as np
import copy
from os import path as path, listdir
import sys
import emcee

from matplotlib.pyplot import close as closeplot
from matplotlib import pyplot as plt

#import corner
import time

def fitEphem(eclFile, T0=None, period=None):
    # data T0, if no T0 is given to us
    if T0 == None:
        print("  No prior T0, using first eclipse in data.")
        T0 = np.min(ts)
    
    print("  Fitting these eclipse times:")
    for t in tl:
        print("  {:.7f}+/-{:.7f} from {}".format(t[0], t[1], t[2]))
    print("\n  Starting from an initial ephem of T0: {}, P: {}".format(T0, period))
    
    def test(params, E):
        # Gets the eclipse number.

        # Extract the params
        T = params[0]
        period = params[1]
        
        # How far are we from the predicted eclipse time
        comp = ((data - T) / period)

        return comp

    def errFunc(p, e, t, t_e):
        diffs = ( test(p, E) - t ) / t_e

        return diffs
    
    out = leastsq(errFunc,
        [T0, period],
        args=(E, ts, t_err),
        full_output=1
    )

    pfinal = out[0]
    covar = out[1]

    P, P_err = pfinal[1], np.sqrt(covar[1][1])
    T0, T0_err = pfinal[0], np.sqrt(covar[0][0])

    print("  Got a T0 of {:.10f}+/-{:.2e}".format(T0, T0_err))
    print("  Got a period of {:.10f}+/-{:.2e}".format(P, P_err))
