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
import glob


class PlotPoints:
    def __init__(self, fig):
        self.fig = fig
        self.xcoords = np.array([])
        self.ycoords = np.array([])
        self.flag = False

    def connect(self):
        self.cidpress = self.fig.canvas.mpl_connect('key_press_event', self.on_press)
        self.cidclick = self.fig.canvas.mpl_connect('button_release_event', self.on_click)
        print("  Hit 'q' to skip these data.\n  Click another button, or the mouse, on initial guesses for ingress and egress:")

    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cidpress)
        self.fig.canvas.mpl_disconnect(self.cidclick)

    def on_press(self, event):
        if 'q' in event.key:
            self.flag = True
            self.disconnect()
            closeplot()
            print("  ")
            return
        print('  added point at %.1f, %.1f' % (event.xdata, event.ydata))
        self.xcoords = np.append(self.xcoords, event.xdata)
        self.ycoords = np.append(self.ycoords, event.ydata)
        if self.xcoords.size == 2:
            self.disconnect()
            closeplot()
            print("  ")

    def on_click(self, event):
        print('  added point at %.1f, %.1f' % (event.xdata, event.ydata))
        self.xcoords = np.append(self.xcoords, event.xdata)
        self.ycoords = np.append(self.ycoords, event.ydata)
        if self.xcoords.size == 2:
            self.disconnect()
            closeplot()
            print("  ")

    def gaussPars(self):
        sep = np.fabs(self.xcoords[0] - self.xcoords[1])
        peak = np.fabs(self.ycoords).mean()
        t0 = self.xcoords.mean()
        sigma = sep/20
        return dict(t0=t0, peak=peak, sep=sep, log_sigma2=np.log(sigma**2))


class TwoGaussians(Model):
    parameter_names = ('t0', 'sep', 'peak', 'log_sigma2')

    def get_value(self, x):
        mu1, mu2 = self.t0 - self.sep/2, self.t0 + self.sep/2
        g1 = np.exp(-0.5 * (x-mu1)**2 * np.exp(-self.log_sigma2))
        g2 = np.exp(-0.5 * (x-mu2)**2 * np.exp(-self.log_sigma2))
        return -self.peak * g1 + self.peak * g2

    def compute_gradient(self, x):
        t, s, p, ls2 = self.t0, self.sep, self.peak, self.log_sigma2

        # some simplifying contractions
        x0 = s/2
        x1 = x - t
        x2 = 0.5*np.exp(-ls2)
        x3 = x2*(x1-x0)**2
        x4 = np.exp(-x3)
        x5 = x0 + x1
        x6 = x2*(x - t + s/2)**2
        x7 = np.exp(-x6)
        x8 = p*x2
        x9 = x7*x8
        x10 = x4*x8
        x11 = 2*(t-x)

        # derivatives, from SymPy
        dp = np.exp(-x3) - np.exp(-x6)
        ds = -x10 * (t - x + x0) + x5*x9
        dt = -x10 * (s + x11) + x9*(x11 - s)
        dls2 = p*x3*x4 - p*x6*x7
        return np.array([dt, ds, dp, dls2])

def tcorrect(tseries, star, observatory, type='B'):
    """
    Correct for light travel time.

    Arguments:
    ----------
    tseries: hipercam.hlog.Tseries
        Time series object

    star: astropy.coordinate.SkyCoord
        Location of star on Sky

    observatory: string
        Observatory name. See coord.EarthLocation.get_site_names() for list

    type: string (default=B)
        Heliocentric (H) or Barcentric (B)

    Returns
    -------
    tseries_corr : hipercam.hlog.Tseries
        Time series object with corrected time axis
    """
    ts = copy.deepcopy(tseries)
    times = Time(tseries.t, format='mjd', scale='utc',
                 location=coord.EarthLocation.of_site(observatory))
    if type == 'B':
        corr = times.light_travel_time(star)
        corr = times.tdb + corr
    else:
        corr = times.light_travel_time(star, 'heliocentric')
        corr = times.utc + corr
    ts.t = corr.mjd
    return ts

def smooth_derivative(tseries, med_half_width, box_half_width):
    """
    Calculate a smoothed version of the lightcurve derivative

    First smooth lightcurve wiht a median filter, then smooth
    numerical derivative with boxcar convolution.

    Parameters
    -----------
    tseries: hipercam.hlog.Tseries
        Time series object

    med_half_width: int
        Half-width of median filter

    box_half_width: int
        Half-width of boxcar filter

    Returns
    --------
    x, dy: np.ndarray
        Locations and values of smoothed derivative
    """
    x = tseries.t.copy()
    y = tseries.y.copy()
    yf = medfilt(y, 2 * med_half_width + 1)
    deriv = (yf[1:] - yf[:-1]) / (x[1:] - x[:-1])
    locs = 0.5 * (x[1:] + x[:-1])
    kernel = Box1DKernel(2 * box_half_width + 1)
    return locs, convolve(deriv, kernel)

def get_tseries(logfile, ccdnam, ap_targ, ap_comp):
    log = Hlog.from_ulog(logfile)
    return log.tseries(ccdnam, ap_targ) / log.tseries(ccdnam, ap_comp)

# Define a cost function for MCMC
def log_like(params, y, gp):
    # print(params)
    gp.set_parameter_vector(params)
    return gp.log_likelihood(y)

# Define a cost function for scipy.minimize
def neg_log_like(params, y, gp):
    gp.set_parameter_vector(params)
    return -gp.log_likelihood(y)
def grad_neg_log_like(params, y, gp):
    gp.set_parameter_vector(params)
    return -gp.grad_log_likelihood(y)[1]



def getEclipseTimes(coords, obsname, myLoc=None):
    '''
coords  - "ra dec" - string, needs to be in a format that astropy can interpret.
    ra  - Target Right Ascension, in hours
    dec - Target Declination, in degrees
obsname - Observing location. Currently must be the /name/ of the observatory.
myLoc   - Working directory.


Searches the current directory for a file containing eclipse times, and fits an ephemeris (T0 and period) to it. 
If <analyse_new> is True, it also seraches <myLoc> for log files, and fits the for an eclipse time. 

Searched <myLoc> for a file containing eclipse times. Searches <myLoc> for logfiles, and fits these for their 
eclipse times.

The technique for this is to make a smoothed plot of the numerical gradient, and look for two mirrored peaks - one 
where the lightcurve enters eclipse (showing as a trough in gradient), and one for egress (showing as a peak in 
gradient). Ideally, they will be mirrors of each other, with the same width and height (though one will be the negative
of the other). 

A double gaussian is fitted to it using a gaussian process, and the midpoint between their peaks is taken to be the 
eclipse time. To characterise the error of the eclipse time, an MCMC is used to sample the found fit. This is beefy, 
and takes a while, but the Hessian we were getting out of scipy.optimize was heavily dependant on initial conditions,
so was untrustworthy.
'''
    ### VARIABLES ###
    ### ------------------------------------------------- ###

    star = coord.SkyCoord(
        coords,
        unit=(u.hour, u.deg)
    )

    if myLoc == None:
        print("  Defaulting to current directory")
        myLoc = path.curdir

    oname = 'eclipse_times.txt'
    oname = '/'.join([myLoc, oname])
    if not path.isfile(oname):
        print("  Couldn't find previous eclipse times file, '{}'. Creating that file.".format(oname))
    
    ### ------------------------------------------------- ###
    
    # tecl list
    tl = []
    if path.isfile(oname):
        print("  Found prior eclipses in '{}'. Using these in my fit.".format(oname))
        # Read in the eclipsetimes.txt file
        fname = '/'.join([myLoc, 'eclipse_times.txt'])
        source_key = {}
        tl = []
        with open(fname, 'r') as f:
            for line in f:
                if line[0] == '#':
                    line = line[1:].strip().split(",")
                    source_key[line[1]] = line[0]
                else:
                    line = line.split(',')
                    line = [float(x) for x in line]
                    line[3] = source_key[line[3]]
                    tl.append()
        for t in tl:
            print("cycle {} -- {:.7f}+/-{:.7f} from {}".format(t[0], t[1], t[2], source_key[t[3]]))


    print("  Grabbing log files...")
    fnames = list(glob.iglob('{}/**/*.log'.format(myLoc), recursive=True))
    
    if len(fnames) == 0:
            print("  I couldn't find any log files! For reference, I searched the following:")
            print("   - {}".format(myLoc))
            exit()
    # List the files we found
    print("  Found these log files: ")
    for i, fname in enumerate(fnames):
        print("  {:2d} - {}".format(i, fname))
    print('  ')
    
    for lf in fnames:
        # lets make the file reading more robust
        try:
            log = Hlog.from_ascii(lf)
        except Exception:
            log = Hlog.from_ulog(lf)
            
        # Get the g band lightcurve, and correct it to the barycentric time
        gband = log.tseries('2', '1') / log.tseries('2', '2')
        gband_corr = tcorrect(gband, star, obsname)
        # Discard the first 10 observations, as they're often junk
        gband_corr = gband_corr[10:]

        x, y = smooth_derivative(gband_corr, 9, 5)
        yerr = 0.001*np.ones_like(x)

        fig, ax = plt.subplots()
        plt.plot(x, y)
        gauss = PlotPoints(fig)
        gauss.connect()
        plt.show()

        if gauss.flag:
            print("  No eclipse taken from these data.")
            continue

        kwargs = gauss.gaussPars()
        # hold values close to initial guesses
        bounds = dict(
            t0=(kwargs['t0']-kwargs['sep']/8, kwargs['t0']+kwargs['sep']/8),
            sep=(0.9*kwargs['sep'], 1.1*kwargs['sep']),
            log_sigma2=(np.log(kwargs['sep']**2/10000), np.log(kwargs['sep']**2/25)),
            peak=(0.9*kwargs['peak'], 1.1*kwargs['peak'])
        )
        kwargs['bounds'] = bounds

        mean_model = TwoGaussians(**kwargs)

        mean, median, std = sigma_clipped_stats(y)
        delta_t = np.mean(np.diff(x))*5
        kernel = terms.RealTerm(log_a=np.log(std**2), log_c=-np.log(delta_t))
        gp = celerite.GP(kernel, mean=mean_model, fit_mean=True)
        gp.compute(x, yerr)
        # print("  Initial log-likelihood: {0}".format(gp.log_likelihood(y)))


        # Fit for the maximum likelihood parameters
        initial_params = gp.get_parameter_vector()
        bounds = gp.get_parameter_bounds()
        

        # Find a solution using Stu's method
        soln = minimize(neg_log_like, initial_params, jac=grad_neg_log_like,
                        method="L-BFGS-B", bounds=bounds, args=(y, gp))
        if not soln.success:
            print('  Warning: may not have converged')
            print(soln.message)

        gp.set_parameter_vector(soln.x)
        mean_model.set_parameter_vector(gp.get_parameter_vector()[2:])
        
        out = soln['x']
        t_ecl = out[2]

        print("  Using MCMC to characterise error at peak likelihood...")


        # Use an MCMC model, starting from the solution we found, to model the errors
        ndim     = 6
        nwalkers = 100

        # Initial positions. Scatter by 0.00001, as this is one above the order of magnitude of the error
        #  we expect on t_ecl.
        p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
        scatter = 0.0001/t_ecl
        p0 *= scatter
        p0 += 1. - scatter
        p0 = np.transpose(np.repeat(out, nwalkers).reshape((ndim, nwalkers))) *p0
        
        # Construct a sampler
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_like, args=[y, gp], threads=1)


        width=40

        # Burn in
        print("")
        nsteps = 100
        for i, result in enumerate(sampler.sample(p0, iterations=nsteps)):
            n = int((width+1) * float(i) / nsteps)
            sys.stdout.write("\r  Burning in...    [{}{}]".format('#'*n, ' '*(width - n)))
        pos, prob, state = result
        
        # Data
        sampler.reset()
        nsteps = 300

        for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
            n = int((width+1) * float(i) / nsteps)
            sys.stdout.write("\r  Sampling data... [{}{}]".format('#'*n, ' '*(width - n))) 
        print("")

        # corner.corner(sampler.flatchain, labels=['???', '???', 't_ecl', '???', 'a', 'b'])
        # plt.show()

        t_ecl = np.mean(sampler.flatchain[:,2])
        err = np.std(sampler.flatchain[:,2])
        sep = np.mean(sampler.flatchain[:,3])

        print("    Got a solution: {:.7f}+/-{:.7f}\n".format(t_ecl, err))

        # print("  Got a Jacobian,\n {}".format(soln['jac']))
        # print("  Got a Hessian,\n {}".format(soln['hess_inv'].todense()))
        # print("  Final log-liklihood: {}".format(soln.fun))

        locflag = input("  What is the source of these data: ")
        tl.append([float(t_ecl), float(err), locflag])

        # Make the maximum likelihood prediction
        mu, var = gp.predict(y, x, return_var=True)
        std = np.sqrt(var)

        # Plot the data
        color = "#ff7f0e"
        fig, ax = plt.subplots(2, 1, sharex=True)
        ax[0].plot(x, y, '.')
        ax[0].plot(x, mu, color=color)
        ax[0].fill_between(x, mu+std, mu-std, color=color, alpha=0.3, edgecolor="none")
        ax[0].plot(x, mean_model.get_value(x), 'k-')
        ax[0].axvline(t_ecl, color='magenta')
        ax[0].set_xlim(left=t_ecl-(1*sep), right=t_ecl+(1*sep))
        
        gband_corr.mplot(ax[1])

        ax[0].set_title("maximum likelihood prediction - {}".format(lf.split('/')[-1]))
        plt.tight_layout()
        plt.show()
    print("  \nDone all the files!")

    # make a key for the data sources
    key = ''
    key_dict = {}
    i = 0
    for t, t_e, source in tl:
        print("source: {}".format(source))
        if source not in key:
            key += "#{},{}\n".format(source, i)
            key_dict[source] = i
            i += 1

    # Sort the list
    tl = sorted(tl,key=lambda x: (x[0]))

    with open(oname, 'w') as f:
        f.write(key)
        for t, t_e, source in tl:
            f.write("<ECLIPSE NUMBER>, {}, {},{}\n".format(t, t_e, key_dict[source]))
    print("  Wrote eclipse data to {}\n".format(oname))

    #TODO:
    # Temporary placeholder. Think about this.
    # - Get the rounded ephemeris fit from the period and T0 supplied?
    # - Might be best to force the user to do this manually, to make it more reliable?
    print("  Please open the file, and edit in the eclipse numbers for each one.\n  Hit enter when you've done this!")
    input()
