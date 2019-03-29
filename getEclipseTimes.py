from hipercam.hlog import Hlog

from astropy import coordinates as coord, units
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
import glob

import mcmc_utils
try:
    from logger import printer
except:
    def printer(string, end='\n'):
        print(string, end='\n')

class PlotPoints:
    def __init__(self, fig):
        self.fig = fig
        ax = self.fig.get_axes()
        self.grad = ax[0]
        self.data = ax[1]
        self.xcoords = np.array([])
        self.ycoords = np.array([])
        self.flag = False

    def connect(self):
        self.cidpress = self.fig.canvas.mpl_connect('key_press_event', self.on_press)
        print("\n\n  Hit 'q' to skip these data.\n  Hit 'a' on initial guesses for ingress and egress\n  +/- set upper and lower limits for the fit\n")

    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cidpress)

    def on_press(self, event):
        if 'q' in event.key:
            self.flag = True
            self.disconnect()
            closeplot()
            print("  ")
            return

        if '-' in event.key:
            print("  lower limit of data is {:.6f}".format(event.xdata))
            self.lowerlim = event.xdata
            self.grad.axvline(self.lowerlim, color='red', linestyle='--')
            self.data.axvline(self.lowerlim, color='red', linestyle='--')

        if '+' in event.key:
            print("  upper limit of data is {:.6f}".format(event.xdata))
            self.upperlim = event.xdata
            self.grad.axvline(self.upperlim, color='green', linestyle='--')
            self.data.axvline(self.upperlim, color='green', linestyle='--')

        if 'a' in event.key:
            print('  added point at {:.1f}, {:.1f}'.format(event.xdata, event.ydata))
            self.xcoords = np.append(self.xcoords, event.xdata)
            self.ycoords = np.append(self.ycoords, event.ydata)
            self.grad.scatter(event.xdata, event.ydata, marker='x', color='black')

        if 'r' in event.key:
            self.xcoords = np.array([])
            self.ycoords = np.array([])

        if self.xcoords.size == 2:
            self.disconnect()
            closeplot()
            print("  ")

    def gaussPars(self):
        sep = np.fabs(self.xcoords[0] - self.xcoords[1])
        peak = np.fabs(self.ycoords).mean()
        t0 = self.xcoords.mean()
        sigma = sep/20
        printer("Took an initial guess with the following parameters:")
        printer("        T0: {}".format(t0))
        printer("      peak: {}".format(peak))
        printer("       sep: {}".format(sep))
        printer("log_sigma2: {}".format(np.log(sigma**2)))
        printer("")
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
        Observatory name. See coord.EarthLocation.get_site_names() for list. If not in the list, assumed
        to be "lat, lon", comma separated.

    type: string (default=B)
        Heliocentric (H) or Barcentric (B)

    Returns
    -------
    tseries_corr : hipercam.hlog.Tseries
        Time series object with corrected time axis
    """
    ts = copy.deepcopy(tseries)
    try:
        location = coord.EarthLocation.of_site(observatory)
    except:
        lat, lon = observatory.split(',')
        print("Attempting to get the earth location from latitude and longitude")
        location = coord.EarthLocation.from_geodetic(lat=lat, lon=lon)

    times = Time(tseries.t, format='mjd', scale='utc', location=location)

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

def read_ecl_file(fname):
    '''
    Reads my eclipse files. Returns an array where each row (tl[i, :]) containing the data for a single eclipse.
    Also returns a key for decoding the source ID


    Arguments:
    ----------
    fname: string
        The eclipse file to be read in

    Returns:
    --------
    source_key: dict
        A dictionary, where each key ( str(int) ) corresponds to the source (string)

    tl: list
        A 2D list containing the cycle number, eclipse time, their error, and source key for each data
    '''
    # tecl list
    tl = []
    if path.isfile(fname):
        # Read in the eclipsetimes.txt file
        source_key = {}
        tl = []
        nObs = 0
        with open(fname, 'r') as f:
            for line in f:
                line = line.strip()
                if line == '':
                    pass
                elif line[0] == '#':
                    line = line[1:].strip().split(",")
                    source_key[line[1]] = line[0]
                else:
                    line = line.split(',')
                    line[:3] = [float(x) for x in line[:3]]
                    line[0] = int(line[0])
                    line[3]  = int(line[3])
                    tl.append(line)
                    nObs += 1
        printer("Found {} prior eclipse times:".format(nObs))
        for e, t, t_err, source in tl:
            printer("-> Cycle: {:5d} -- {:.7f}+/-{:.7f} from {}".format(e, t, t_err, source_key[str(source)]))
    else:
        printer("Could not find the file '{}' to read eclipse data from! Creating a new file.".format(fname))
        source_key, tl = {}, []
    printer("")
    return source_key, tl

def write_ecl_file(source_key, tl, oname):
    """
    Saves eclipse data to a file.


    Arguments:
    ----------
    source_key: dict
        A dict, with keys of '1', '2', '3', etc, where each key corresponds to a source string

    tl: numpy.array
        An array of shape (N, 4), containing the cycle number, eclipse time, time error, and source ID
        of N eclipses

    Returns:
    --------
    None
    """
    # make a key for the data sources
    key = ''
    i = 0
    for c, t, t_e, source in tl:
        source = source_key[str(source)]
        if source not in key:
            key += "#{},{}\n".format(source, i)
            i += 1

    # Sort the list
    tl = sorted(tl,key=lambda x: (x[1]))

    with open(oname, 'w') as f:
        f.write(key)
        for c, t, t_e, source in tl:
            f.write("\n{},{},{},{}".format(c, t, t_e, source))
    printer("Wrote eclipse data to {}\n".format(oname))

    return


def getEclipseTimes(coords, obsname, myLoc=None):
    '''
    Searches <myLoc> for .log files, and uses them to get the times of the eclipses.

    The technique for this is to make a smoothed plot of the numerical gradient, and look for two mirrored peaks - one
    where the lightcurve enters eclipse (showing as a trough in gradient), and one for egress (showing as a peak in
    gradient). Ideally, they will be mirrors of each other, with the same width and height (though one will be the negative
    of the other).

    A double gaussian is fitted to it using a gaussian process, and the midpoint between their peaks is taken to be the
    eclipse time. To characterise the error of the eclipse time, an MCMC is used to sample the found fit. This is beefy,
    and takes a while, but the Hessian error matrix we were getting out of scipy.optimize was heavily dependant on initial
    conditions, so was untrustworthy.


    Arguments:
    ----------
    coords: str
        The RA and Dec of the stars in the eclipses you're fitting. Note that all data being fitted is assumed to be for
        the same object, hence the RA and Dec used in each log file is the same. i.e., make sure you're not fitting data
        for more than one object at once!
        Note: must be readable by astropy!

    obsname: str
        The observatory name. See coord.EarthLocation.get_site_names() for a list. If a site is not in the registry,
        this string is assumed to be longitude and latitude, and will be attempted again.

    myLoc: str, default None
        The directory to search for eclipses. If None, searches the current working directory.

    Returns:
    --------
    None, but creates a file with eclipse times in it.
    '''
    plt.ion()
    printer("\n\n--- Getting eclipse times from the data ---")

    star = coord.SkyCoord(
        coords,
        unit=(units.hour, units.deg)
    )

    # Where are we working?
    if myLoc == None:
        myLoc = path.curdir
        printer("Defaulting to current directory: {}".format(myLoc))

    # Where am I looking for prior data, and saving my new data?
    oname = 'eclipse_times.txt'
    oname = '/'.join([myLoc, oname])
    source_key, tl = read_ecl_file(oname)

    # What am I using to get new data from?
    printer("Grabbing log files...")
    fnames = list(glob.iglob('{}/**/*.log'.format(myLoc), recursive=True))
    fnames = sorted(fnames)

    if len(fnames) == 0:
        printer("I couldn't find any log files in:")
        printer("{}".format(myLoc))
        raise FileNotFoundError

    # List the files we found
    printer("Found these log files: ")
    for i, fname in enumerate(fnames):
        printer("  {:>2d} - {}".format(i, fname))
    printer('  ')

    for lf in fnames:
        # lets make the file reading more robust

        try:
            log = Hlog.from_ascii(lf)
            if log == {}:
                raise Exception
        except Exception:
            printer("Using the ulog funtion to read data...")
            log = Hlog.from_ulog(lf)
        aps = log.apnames

        printer("File: {}".format(lf))
        if len(aps['1']) < 2:
            printer("-> Not looking for eclipses in {}, as only one aperture in the file.".format(lf))
            continue

        # Get the first CCD lightcurve, and correct it to the barycentric time
        inspect = log.tseries('1', '1') / log.tseries('1', aps['1'][1])
        inspect_corr = tcorrect(inspect, star, obsname)
        # Discard the first 10 observations, as they're often junk
        inspect_corr = inspect_corr[10:]

        x, y = smooth_derivative(inspect_corr, 9, 5)
        yerr = 0.001*np.ones_like(x)

        fig, ax = plt.subplots(2, figsize=[16,8])
        ax[0].set_title("{}".format(lf))
        ax[0].plot(x, y)

        ax[1].set_title('Lightcurve:')
        inspect_corr.mplot(ax[1])

        gauss = PlotPoints(fig)
        gauss.connect()
        plt.tight_layout()
        plt.show(block=True)

        try:
            lowerlim = gauss.lowerlim
        except:
            lowerlim = x.min()

        try:
            upperlim = gauss.upperlim
        except:
            upperlim = x.max()

        # Apply upper/lower limits
        mask = (x < upperlim) * (x > lowerlim)
        mask = np.where(mask==1)

        y    = y[mask]
        yerr = yerr[mask]
        x    = x[mask]

        if gauss.flag:
            printer("-> No eclipse taken from {}".format(lf))
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


        # Find a solution using Stu's minimisation method
        soln = minimize(neg_log_like, initial_params, jac=grad_neg_log_like,
                        method="L-BFGS-B", bounds=bounds, args=(y, gp))
        if not soln.success:
            printer('  Warning: may not have converged')
            printer(soln.message)

        gp.set_parameter_vector(soln.x)
        mean_model.set_parameter_vector(gp.get_parameter_vector()[2:])

        out = soln['x']
        t_ecl = out[2]

        printer("Using MCMC to characterise error at peak likelihood...")


        # Use an MCMC model, starting from the solution we found, to model the errors
        ndim     = 6
        nwalkers = 100

        # Initial positions.
        p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
        scatter = 0.0005 # Scatter by ~ 40s in time
        p0 *= scatter
        p0 -= (scatter/2)
        p0 = np.transpose(np.repeat(out, nwalkers).reshape((ndim, nwalkers))) + p0

        # print(86400*(p0[:,2] - np.transpose(np.repeat(out, nwalkers).reshape((ndim, nwalkers)))[:,2]))

        # Construct a sampler
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_like, args=[y, gp], threads=1)


        width=40

        # Burn in
        print("")
        nsteps = 1000
        try:
            for i, result in enumerate(sampler.sample(p0, iterations=nsteps)):
                n = int((width+1) * float(i) / nsteps)
                # print(result[0])
                sys.stdout.write("\r  Burning in...    [{}{}]".format('#'*n, ' '*(width - n)))
            pos, prob, state = result

            # Data
            sampler.reset()
            nsteps = 2000

            for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
                n = int((width+1) * float(i) / nsteps)
                sys.stdout.write("\r  Sampling data... [{}{}]".format('#'*n, ' '*(width - n)))
            print("")

            # chain = sampler.flatchain
            # fig = mcmc_utils.thumbPlot(chain,['g1', 'g2', 'T0', 'sep', 'peak', 'log_sigma2'])
            # fig.savefig('/'.join([myLoc, 'eclipse_{}_cornerPlot.pdf'.format(lf.split('/')[-1])]))
            # plt.show(block=False)

            t_ecl = np.mean(sampler.flatchain[:,2])
            err = np.std(sampler.flatchain[:,2])
            sep = np.mean(sampler.flatchain[:,3])

            printer("Got a solution: {:.7f}+/-{:.7f}\n".format(t_ecl, err))

            # Make the maximum likelihood prediction
            mu, var = gp.predict(y, x, return_var=True)
            std = np.sqrt(var)

            # Plot the data
            color = "#ff7f0e"
            plt.close('all')
            fig, ax = plt.subplots(2, 1, sharex=True)
            ax[0].plot(x, y, '.')
            ax[0].plot(x, mu, color=color)
            ax[0].fill_between(x, mu+std, mu-std, color=color, alpha=0.3, edgecolor="none")
            ax[0].plot(x, mean_model.get_value(x), 'k-')
            ax[0].axvline(t_ecl, color='magenta')

            inspect_corr.mplot(ax[1])
            ax[1].set_title('Lightcurve')
            ax[1].axvline(t_ecl, color='magenta')
            ax[1].axvline(t_ecl+(sep/2.), color='red')
            ax[1].axvline(t_ecl-(sep/2.), color='red')

            ax[0].set_title("maximum likelihood prediction - {}".format(lf.split('/')[-1]))
            plt.tight_layout()
            print("  Plotting fit...")
            plt.show(block=False)

            cont = input("  Save these data? y/n: ")
            if cont.lower() == 'y':
                locflag = input("    What is the source of these data: ")

                key = '-1' # This ensures that if source_key is empty, the new data are pushed to index '0'
                for key in source_key:
                    if locflag == source_key[key]:
                        locflag = key
                        break
                if locflag != key:
                    key = str(int(key)+1)
                    source_key[key] = locflag
                    locflag = key
                tl.append(['0', float(t_ecl), float(err), locflag])
                printer("Saved the data: {}".format(['0', float(t_ecl), float(err), locflag]))
            else:
                printer("  Did not store eclipse time from {}.".format(lf))
            plt.close()
            printer("")
        except celerite.solver.LinAlgError:
            printer('  Celerite failed to factorize or solve matrix. This can happen when the data are poorly fitted by the double gaussian!')
            printer("  Skipping this file.")

    printer("\nDone all the files!")

    write_ecl_file(source_key, tl, oname)
    plt.ioff()

    #TODO:
    # Temporary placeholder. Think about this.
    # - Get the rounded ephemeris fit from the period and T0 supplied?
    # - Might be best to force the user to do this manually, to make it more reliable?
    printer("This string might help:\ncode {}".format(path.abspath(oname)))
    printer("Please open the file, and edit in the eclipse numbers for each one.")
    input("Hit enter when you've done this!")
