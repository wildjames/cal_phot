from hipercam.hlog import Hlog

from astropy import coordinates as coord, units as u
from astropy.time import Time
from astropy.convolution import Box1DKernel, convolve
from astropy.stats import sigma_clipped_stats

from scipy.signal import medfilt
from scipy.optimize import minimize, leastsq as lsq

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

def fitEphem(myLoc, T0, period):
    # Read in the eclipsetimes.txt file
    fname = '/'.join([myLoc, 'eclipse_times.txt'])
    source_key = {}
    tl = []
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip()
            if line[0] == '#':
                line = line[1:].strip().split(",")
                source_key[line[1]] = line[0]
            else:
                line = line.split(',')
                line[:3] = [float(x) for x in line[:3]]
                line[0] = int(line[0])
                line[3]  = source_key[line[3]]
                tl.append(line)

    print("  Fitting these eclipse times:")
    for t in tl:
        print("  Cycle: {:5d} -- {:.7f}+/-{:.7f} from {}".format(t[0], t[1], t[2], t[3]))
    print("\n  Starting from an initial ephem of T0: {}, P: {}".format(T0, period))
    
    ### Fitting

    # def test(params, E):
    #     # Gets the eclipse number.

    #     # Extract the params
    #     T = params[0]
    #     period = params[1]
        
    #     # How far are we from the predicted eclipse time
    #     calc = T + (E*period)

    #     return calc

    # def errFunc(p, data):
    #     e = np.array([i[0] for i in data])
    #     t = np.array([i[1] for i in data])
    #     t_e = np.array([i[2] for i in data])

    #     diffs = ( test(p, e) - t ) / t_e

    #     return diffs
    
    # out = leastsq(errFunc,
    #     [T0, period],
    #     args=(tl),
    #     full_output=1
    # )

    # pfinal = out[0]
    # covar = out[1]

    # T0, T0_err = pfinal[0], np.sqrt(covar[0][0])
    # P, P_err   = pfinal[1], np.sqrt(covar[1][1])


    def fitfunc(p,x):
        return p[0] + p[1]*x
    def errfunc(p,x,y,err):
        return (y-fitfunc(p,x)) / err

    ### params = numpy.array([55202.0579,0.05298466])

    params = [T0, period]

    x =  np.array([i[0] for i in tl])
    y =  np.array([i[1] for i in tl])
    ey = np.array([i[2] for i in tl])

    out = lsq(errfunc,params,args=(x,y,ey),full_output=1)
    pfinal = out[0]
    covar = out[1]

    P, P_err = pfinal[1], np.sqrt(covar[1][1])
    T0, T0_err = pfinal[0], np.sqrt(covar[0][0])

    chisq = (y - fitfunc(pfinal, x))**2 / (ey**2)
    chisq /= 

    ### Reporting
    print("  Got a T0 of {:.10f}+/-{:.2e}".format(T0, T0_err))
    print("  Got a period of {:.10f}+/-{:.2e}".format(P, P_err))
    print("This fit had a chisq value of {}".format(chisq))
    print('')
    print("  (Obs) - (Calc), sec | Cycle Number")
    for t in tl:
        dT = fitfunc(pfinal, t[0]) - t[1]
        dT *= 24*60*60
        print(" {:>20.4f} | {:d}".format(dT , t[0] ))
    


    return T0, P