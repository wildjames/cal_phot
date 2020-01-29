import matplotlib.pyplot as plt
import numpy as np
from astropy import coordinates as coord
from astropy import time
from astropy import units as u
from astropy.coordinates import AltAz
from hipercam.hlog import Hlog
from scipy import optimize as opt


# The model to fit
def model(airmass, k, F0):
    # If F0 or k is negative, we're in bad territory
    if F0 < 0.0 or k < 0.0:
        return -np.inf
    pred = -2.5*np.log10(F0) + k * airmass
    if np.any(np.isnan(pred)):
        pred = -np.inf
    return pred

def chisq_func(model_args, model_function, x_data, y_data, yerr=None, ln_like=False):
    '''Call the model function on the x_data, compute chisq on the observations.
    Model takes (x, *model_args)'''

    # I need model_args to be a toople
    model_args = tuple(model_args)
    model_data = model_function(x_data, *model_args)

    if yerr is None:
        yerr = np.ones_like(y_data)

    chisq = np.sum( ((y_data - model_data) / yerr)**2 )

    # emcee wants to maximise the diagnostic it's given.
    #   maximizing the negative chisq -> minimizing the normal chisq.
    # This assumes that the noise on the data is gaussian. Otherwise, technically invalid.
    if ln_like:
        return -0.5*chisq

    return chisq

def scipy_opt(model_func, initial_guess, data, *args, **kwargs):
    '''Scipy convenience function. *args and **kwargs are passed to scipy.optimise.'''
    # Construct the argument tuple for chisq_func
    chisq_args = (model_func, *data)

    # Fit with scipy
    result = opt.minimize(chisq_func, initial_guess, args=chisq_args, *args, **kwargs)
    return result

def fit_airmass(target_loc, observatory, data, CCD, aperture, quiet=True):
    '''
    Fit extinction coefficients to an observation

    Inputs
    ------
      - target_loc: astropy.coordinates.SkyCoord
        - The location of the target star in the sky
      - observatory: astropy.coordinates.EarthLocation
        - The location of the observing site
      - data: hipercam.hlog.Hlog
        - Target data to be fitted. Will have outliers removed
      - CCD: str
        - The CCD to fit
      - aperture: str
        - The aperture in that CCD to use
      - quiet: bool
        - Whether to report more stuff or not

    Output
    ------
      - fit: scipy.optimize.OptimizeResult
        - The result of the fit
    '''
    if not quiet:
        print("The data has apertures:")
        print("  {}".format(data.apnames))
        print("Using CCD {}, aperture {}".format(CCD, aperture))

    # Grab the data I care about
    CCD = str(CCD)
    aperture = str(aperture)
    star = data.tseries(CCD, aperture)

    # Filter out nans an 0.0 data
    star = star.remove_outliers()

    # TODO: Filter out cloudy observations

    # Covert to electons/second, then to magnitudes
    exp_time = np.mean(data[CCD]['Exptim'])
    electronic_mags = -2.5 * np.log10(star.y/exp_time)
    electronic_mags_err = (np.log(np.e)/np.log(10)) * (star.ye/star.y)
    if not quiet:
        print("\nThe error in magnitude will be {:.3f} * (err/flx)".format(np.log(np.e)/np.log(10)))

    ## Calculate the airmasses of each frame
    # obs_T can be an array
    obs_T = star.t
    obs_T = time.Time(obs_T, format='mjd')
    if not quiet:
        print("\n\nThe observations start at {}".format(obs_T[0].iso))
        print("The observations end at   {}".format(obs_T[-1].iso))

    # Define the altAz frame, containing the time and location
    alt_az_frame = AltAz(
        obstime=obs_T,
        location=observatory
    )

    # Transform the star location at that time
    target_loc_AltAz = target_loc.transform_to(alt_az_frame)
    if not quiet:
        print("\n\nTarget altitude starts at: {:.3f}".format(target_loc_AltAz[0].alt))
        print("Target altitude ends at:   {:.3f}".format(target_loc_AltAz[-1].alt))

    # Get the object airmasses
    zenith_angle = 90. - target_loc_AltAz.alt.deg
    zenith_angle_rad = np.deg2rad(zenith_angle)
    airmasses = 1. / np.cos(zenith_angle_rad)

    if not quiet:
        print("\n\nInitial airmass: {:.3f}".format(airmasses[0]))
        print("Final airmass: {:.3f}".format(airmasses[-1]))
        # Plot the graphs
        fig, ax = plt.subplots(figsize=(12,8))
        ax.set_title("Raw electron counts")

        times_hours = (star.t - star.t[0]) * 60*24
        ax.errorbar(
            times_hours, star.y, yerr=star.ye,
            fmt='x ', color='black'
        )

        ax.set_ylabel("Electron counts, C_{obs}")
        ax.set_xlabel("Time, minutes")
        plt.show()

        fig, ax = plt.subplots(figsize=(12,8))
        ax.set_title("Converted to magnitudes and airmasses")
        ax.errorbar(
            airmasses, electronic_mags, yerr=electronic_mags_err,
            fmt='x ', color='black'
        )
        ax.set_ylabel("$-2.5 log_{10}(C_{obs})$")
        ax.set_xlabel("Airmass, sec(z)")
        plt.show()

    # Fitting
    initial_guess = (
        0.20,  # k
        np.mean(star.y),   # F0
    )

    # Optimize fit
    data = (airmasses, electronic_mags, electronic_mags_err)
    fit = scipy_opt(model, initial_guess, data)
    mean = fit['x']
    err = fit['hess_inv']

    if not quiet:
        print("CCD: {}, aperture: {} gave a value of k: {:.3f}".format(CCD, aperture, mean[0]))
        print("Hessian inv matrix: \n{}".format(err))

    model_airmass = np.linspace(np.min(airmasses), np.max(airmasses), 100)
    model_mags = model(model_airmass, *fit['x'])

    # Plot the fit
    fig, ax = plt.subplots(figsize=(12,8))
    ax.set_title("CCD {}, Ap {} Magnitudes and airmasses, Fitted".format(CCD, aperture))

    ax.errorbar(
        airmasses, electronic_mags, yerr=electronic_mags_err,
        fmt='x ', color='black', label='Data', zorder=2
    )

    ax.plot(
        model_airmass, model_mags, color='red',
        label="k: {:.3f} +/- {:.3f}".format(mean[0], err[0,0]),
        zorder=3
    )

    ax.set_ylabel("$-2.5 log_{10}(C_{obs})$")
    ax.set_xlabel("Airmass, sec(z)")
    ax.legend()
    plt.show()

    return fit

if __name__ == "__main__":
    import argparse
    import pandas as pd

    desc = '''
Takes a log file and calculates an atmospheric extinction coefficient in mags/airmass for each aperture in the file.

There are a few caveats to bear in mind.
  - You need to use an observation on a CLEAR, PHOTOMETRIC night. Otherwise, this will be incorrect! Clouds would contribute to the extinction and give you a wrong value.
  - The larger the airmass range, the better. The logs on `deneb` list the airmass range for an observation, so use that to inform what  run you use.
  - Variable targets will result in junk. There's a lazy way and a smart way to deal with this - use many apertures to make sure that at least some are constant sources, or actually check each target's RA and Dec in catalogues for variability.

The RA and Dec are used to calculate the airmass at whatever time, so the tel RA, Dec from the logs will be fine
'''

    parser = argparse.ArgumentParser(
        description=desc,
        prefix_chars='-'
    )

    parser.add_argument(
        'logfile',
        help='hipercam pipeline log file to use. All apertures will be used to calculate extinction'
    )
    parser.add_argument(
        'RA',
        help='RA of the target frame, hourangle'
    )
    parser.add_argument(
        'Dec',
        help='Declination of the target frame, deg'
    )
    parser.add_argument(
        'observatory',
        help="The observing location. Must be parsable by astropy's of_site()"
    )
    parser.add_argument(
        "--loud",
        help='Print debugging information',
        action="store_true"
    )

    args = parser.parse_args()
    data = args.logfile
    RA, Dec = args.RA, args.Dec
    observatory = args.observatory
    debug = args.loud

    print('log file = {}'.format(data))

    data = Hlog.read(data)

    target_loc = coord.SkyCoord(RA, Dec, unit=(u.hourangle, u.deg))
    print("Target location: {}".format(target_loc))

    obs_loc = coord.EarthLocation.of_site(observatory)
    print("Observing location: {}".format(obs_loc))

    extinction_coefficients = {}
    for CCD, apertures in data.apnames.items():
        print("CCD {}".format(CCD))

        extinction_coefficients[CCD] = []
        for aperture in apertures:
            fit = fit_airmass(target_loc, obs_loc, data, CCD, aperture, quiet=(not debug))
            ext = fit['x'][0]
            extinction_coefficients[CCD].append(ext)
            print("  Aperture {}:".format(aperture))
            print("    Extinction/Airmass = {:.3f}".format(ext))

    print("\n\n\nMeans:")
    for CCD, ks in extinction_coefficients.items():
        print("CCD {}: k = {:.4f}".format(CCD, np.mean(ks)))
