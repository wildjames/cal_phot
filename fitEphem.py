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

import pylab
import mcmc_utils as mu

from getEclipseTimes import read_ecl_file
from logger import printer

#import corner
import time

def fitEphem(myLoc, T0, period, simple=False):
    """
    Takes eclipse time data from a file, and fits ephemeris parameters (T0, period) to them
    from an initial guess.

    There are two fitting options available: simple and complex.
    simple:
        Fits a simple linear regression fit to the data. Report the minimum chisq params
    complex:
        Uses an MCMC chain and a gaussian process to fit the data. Also considers that the
        reported error is not accurate, and allows data from different sources to have their
        error bars scaled to get a chisq closer to unity. This is far more reliable and accurate,
        but takes significantly longer to run.

    Arguments:
    ----------
    myLoc: str
        Tells the function where to look for both logfiles and prior eclipse data

    T0: float
        Initial guess for T0

    period: float
        Initial guess for the period

    simple: bool
        If true, use a linear regression fit on the data. If false, use an MCMC chain.

    Returns:
    --------
    T0: float
        Best T0
    period: float
        Best period
    """
    printer("\n\n--- Fitting ephemeris to data ---")

    # Read in the eclipsetimes.txt file
    fname = '/'.join([myLoc, 'eclipse_times.txt'])
    source_key, tl = read_ecl_file(fname)

    ### Fitting
    x        = np.array([float(x[0]) for x in tl]) # cycle number
    y        = np.array([float(x[1]) for x in tl]) # Time
    ey       = np.array([float(x[2]) for x in tl]) # Error
    obsCodes = np.array([x[3] for x in tl]) # Data source

    params = [T0, period]

    if simple:
        def fitfunc(p,x):
            return p[0] + p[1]*x
        def errfunc(p,x,y,err):
            return (y-fitfunc(p,x)) / err

        out = lsq(errfunc,params,args=(x,y,ey),full_output=1)
        pfinal = out[0]
        covar = out[1]

        P, P_err = pfinal[1], np.sqrt(covar[1][1])
        T0, T0_err = pfinal[0], np.sqrt(covar[0][0])

        chisq = (y - fitfunc(pfinal, x))**2 / (ey**2)
        chisq = sum(chisq) / (len(chisq)-2)

        resy = y - fitfunc(pfinal, x)

    else:
        # Use Stu's version of the fitting, which considers error
        def model(pars,x):
            return pars[0] + pars[1]*x

        def chisq(pars,x,y,yerr):
            resids = ( y - model(pars,x) ) / yerr
            return np.sum(resids*resids)

        def reducedChisq(pars,x,y,yerr):
            return chisq(pars,x,y,yerr) / (len(x) - len(pars) - 1)

        def ln_likelihood(pars,x,y,yerr,obsCodes):
            errs = yerr.copy()
            # scale errors by amount set by observatory code
            emul_factors = pars[2:]
            multipliers  = emul_factors[obsCodes-1]
            errs = errs*multipliers
            return -0.5*(np.sum( np.log( 2.0*np.pi*errs**2 ) ) + chisq(pars,x,y,errs))

        def ln_prior(pars):
            lnp = 0.0
            # only priors are on error scaling - assume good to 2 orders of magnitude
            prior = mu.Prior('log_uniform',0.01,100)
            for param in pars[2:]:
                lnp += prior.ln_prob(param)
            return lnp

        def ln_prob(pars,x,y,yerr,obsCodes):
            lnp = ln_prior(pars)
            if np.isfinite(lnp):
                return lnp + ln_likelihood(pars,x,y,yerr,obsCodes)
            else:
                return lnp

        # p = [T0, period]
        def fitfunc(p,x):
            return p[0] + p[1]*x


        # Get the number of sources, which will each have their errors scaled differently
        numObs = len(source_key)


        # Set up the MCMC chain.
        npars = 2+numObs                # The parameters are the T0, P, and the error scale of each source
        nwalkers = max(16,4*(2+numObs)) # The number of walklers wants to be at least 16,
                                        # or enough to sample our parameter space properly
        params = np.lib.polynomial.polyfit(x,y,1).tolist() # just fit times and cycle numbers to get guesses
        params = params[::-1]

        # Construct a list of parameters for the model, i.e. T0, P, and the error scale factor of each source
        nameList = ['T0','P']
        for i, key in enumerate(source_key):
            # add error factor for each source
            params.append(1.0)
            sourceName = source_key[key].replace(' ', '_')
            nameList.append('obs_err_{:s}'.format(sourceName))
        guessP = np.array(params)

        # Initialise the sampler
        p0 = emcee.utils.sample_ball(guessP,0.01*guessP,size=nwalkers)
        sampler = emcee.EnsembleSampler(nwalkers,npars,ln_prob,args=[x,y,ey,obsCodes],threads=1)

        #burn in
        nburn=1000
        pos, prob, state = mu.run_burnin(sampler, p0,nburn)

        #production
        sampler.reset()
        nprod = 9000
        sampler = mu.run_mcmc_save(sampler,pos,nprod,state,"{}/ephemerisChain.txt".format(myLoc))
        chain = mu.flatchain(sampler.chain,npars,thin=1)

        # Gather and report the best values
        bestPars = []
        for i in range(npars):
            par = chain[:,i]
            lolim,best,uplim = np.percentile(par,[16,50,84])
            printer("{:>20s} = {:.10f} +{:.10f} -{:.10f}".format(nameList[i],best,uplim-best,best-lolim))
            bestPars.append(best)

            if nameList[i] == 'T0':
                T0 = best
                T0_err = uplim-lolim
            if nameList[i] == 'P':
                P = best
                P_err = uplim-lolim
        fig = mu.thumbPlot(chain,nameList)
        fig.savefig('/'.join([myLoc, 'ephemeris_cornerPlot.pdf']))
        printer("Saved a corner plot of the MCMC fit (including error scaling factors) to:\n-> {}".format('/'.join([myLoc, 'ephemeris_cornerPlot.pdf'])))
        plt.close('all')

        resy = 86400.0*(y-model(bestPars,x))
        errs = ey.copy()
        # scale errors by amount set by observatory code
        emul_factors = np.array(bestPars[2:])
        multipliers  = emul_factors[obsCodes-1]
        errs = errs*multipliers
        ey   = 86400.0*errs

        chisq = np.sum(resy**2.0/ey**2.0)
        printer("\nChisq = {:.1f}, with {:d} degrees of freedom".format(chisq, int(x.size - 2)))

    ### Reporting
    printer("Got a T0 of     {:>5.15f}+/-{:<.2e}".format(T0, T0_err))
    printer("Got a period of {:>5.15f}+/-{:<.2e}".format(P, P_err))
    printer("This fit had a reduced chisq value of {:.3f}".format(chisq))
    printer('')
    printer("Source          |  (Obs) - (Calc), sec | Cycle Number")
    for t in tl:
        dT = fitfunc([T0, P], t[0]) - t[1]
        dT *= 24*60*60
        printer(" {:<14s} | {:>20.4f} | {:d}".format(source_key[str(int(t[3]))], dT , int(t[0]) ))

    plt.errorbar(x,resy,yerr=ey,fmt='o',color='k',ecolor='k')
    plt.axhline(ls='--',color='k')
    plt.xlabel('Cycle No.')
    plt.ylabel('O-C (s)')
    plt.show()

    printer("")
    return T0, P