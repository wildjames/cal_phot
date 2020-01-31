from os import path as path

import emcee
import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt
from scipy.optimize import leastsq as lsq

from .mcmc_utils import *
from .getEclipseTimes import read_ecl_file

try:
    from .logger import printer
except ImportError:
    def printer(string, end='\n'):
        print(string, end='\n')


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

    # multipliers  = emul_factors[obsCodes-1]
    multipliers = np.array([1.0 for _ in emul_factors[obsCodes-1]])

    errs = errs*multipliers

    return -0.5*(np.sum( np.log( 2.0*np.pi*errs**2 ) ) + chisq(pars,x,y,errs))

def ln_prior(pars):
    lnp = 0.0
    # only priors are on error scaling - assume good to 2 orders of magnitude
    prior = Prior('log_uniform',0.01,100)
    for param in pars[2:]:
        lnp += prior.ln_prob(param)
    return lnp

def ln_prob(pars,x,y,yerr,obsCodes):
    lnp = ln_prior(pars)
    if np.isfinite(lnp):
        return lnp + ln_likelihood(pars,x,y,yerr,obsCodes)
    else:
        return lnp

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
    fname = path.join(myLoc, 'EPHEMERIS', 'eclipse_times.txt')

    if not path.isfile(fname):
        raise FileNotFoundError("I don't have any eclipse times to fit!")

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
        # Get the number of sources, which will each have their errors scaled differently
        numObs = len(source_key)


        # Set up the MCMC chain.
        npars = 2+numObs                # The parameters are the T0, P, and the error scale of each source
        nwalkers = max(16,4*(2+numObs)) # The number of walklers wants to be at least 16,
                                        # or enough to sample our parameter space properly

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
        sampler = emcee.EnsembleSampler(nwalkers,npars,ln_prob,args=[x,y,ey,obsCodes])

        #burn in
        nburn=5000
        pos, prob, state = run_burnin(sampler, p0,nburn)

        #production
        sampler.reset()
        nprod = 2000

        chain_fname = path.join(myLoc, "EPHEMERIS", "ephemeris_chain.txt")

        sampler = run_mcmc_save(sampler, pos, nprod, state, chain_fname)
        chain = flatchain(sampler.chain, npars, thin=1)

        # Gather and report the best values
        bestPars = []
        for i in range(npars):
            par = chain[:,i]
            lolim,best,uplim = np.percentile(par,[16,50,84])
            # printer("{:>20s} = {:.10f} +{:.10f} -{:.10f}".format(nameList[i],best,uplim-best,best-lolim))
            bestPars.append(best)

            if nameList[i] == 'T0':
                printer("Old T0: {:.8f}".format(T0))
                T0 = best
                T0_err = uplim-lolim
                printer("New T0: {:.8f}+/-{:.8f}\n".format(T0, T0_err))
            if nameList[i] == 'P':
                printer("Old P:  {:.8f}".format(period))
                P = best
                P_err = uplim-lolim
                printer("New P: {:.8f}+/-{:.8f}\n".format(P, P_err))

        fig = thumbPlot(chain,nameList)

        corner_fname = path.join(myLoc, "EPHEMERIS", "ephemeris_corner_plot.pdf")
        fig.savefig(corner_fname)

        printer("Saved a corner plot of the MCMC fit (including error scaling factors) to:\n-> {}".format(corner_fname))
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

    def fitfunc(p,x):
        return p[0] + p[1]*x
    for t in tl:
        dT = fitfunc([T0, P], t[0]) - t[1]
        dT *= 24*60*60
        printer(" {:<14s} | {:>20.4f} | {:d}".format(source_key[str(int(t[3]))], dT , int(t[0]) ))

    # Each error code wants to be a different color
    codes = set(obsCodes) # set() strips duplicates
    codes = list(codes)   # list() allows indexing
    CB_color_cycle = ['black', '#ff7f00', '#4daf4a',
                      '#377eb8', '#a65628', '#984ea3',
                      '#999999', '#e41a1c', '#dede00']
    colors = []
    for code in obsCodes:
        i = codes.index(code)
        colors.append(CB_color_cycle[i])

    plt.errorbar(x,resy,yerr=ey,marker='',ecolor='k', linestyle='', zorder=1)
    plt.scatter(x, resy, c=colors, marker='o', zorder=2)
    plt.axhline(ls='--',color='k', zorder=1)
    plt.xlabel('Cycle No.')
    plt.ylabel('O-C (s)')

    scatter_fname = path.join(myLoc, "EPHEMERIS", "ephemeris_scatter.pdf")
    ephem_fname = path.join(myLoc, "EPHEMERIS", "ephemeris_fit.txt")
    with open(ephem_fname, 'w') as f:
        f.write("Got a T0 of     {:>5.15f}+/-{:<.2e}\n".format(T0, T0_err))
        f.write("Got a period of {:>5.15f}+/-{:<.2e}\n".format(P, P_err))
        f.write("This fit had a reduced chisq value of {:.3f}\n".format(chisq))
        f.write('\n')
        f.write("Source          |  (Obs) - (Calc), sec | Cycle Number\n")
        for t in tl:
            dT = fitfunc([T0, P], t[0]) - t[1]
            dT *= 24*60*60
            f.write(" {:<14s} | {:>20.4f} | {:d}\n".format(source_key[str(int(t[3]))], dT , int(t[0]) ))

    plt.savefig(scatter_fname)
    plt.show()

    printer("")
    return T0, P
