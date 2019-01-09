import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import copy
from astropy import time, coordinates as coord, units as u
from astropy.coordinates import AltAz

import hipercam as hcam
from constructReference import construct_reference, get_comparison_magnitudes
from getEclipseTimes import read_ecl_file
from logger import printer

def sdss_mag2flux(mag):
    '''Takes an SDSS magnitude, returns the corresponding flux in [mJy]'''
    alpha = 3631e3

    flux = 10**(-mag/2.5)
    flux*= alpha

    return flux

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
    times = time.Time(tseries.t, format='mjd', scale='utc',
                 location=coord.EarthLocation.of_site(observatory))
    if type == 'B':
        corr = times.light_travel_time(star)
        corr = times.tdb + corr
    else:
        corr = times.light_travel_time(star, 'heliocentric')
        corr = times.utc + corr
    ts.t = corr.mjd
    return ts

def calc_E(T, T0, P):
    E = (T-T0) / P
    return E
def calc_E_Err(T, T0, P, T_err, T0_err, P_err):
    N = E(T, T0, P)
    err = np.sqrt(
        ( np.sqrt(T_err**2 + T0_err**2)/(T-T0) )**2 +
        (P_err/P)**2
    )
    E_err = N * err
    return E_err


def combineData(oname, coords, obsname, T0, period, SDSS=True, std_fname=None, comp_fnames=None,
                binsize=1, myLoc='.', ext=None, fnames=None, std_coords=None, std_mags=None):
    '''
    Takes a set of *CAM observations (and data on the system), and produces a set of phase folded lightcurves.

    If we're in the SDSS field, each .log file needs a corresponding .coords file that contains the RA and Dec of each aperture:
        <CCD1 ap1 RA> <Dec>
        <CCD1 ap2 RA> <Dec>

        <CCD2 ap1 RA> <Dec>
        <CCD2 ap2 RA> <Dec>

        <CCD3 ap1 RA> <Dec>
        <CCD3 ap2 RA> <Dec>
        <CCD3 ap3 RA> <Dec>

    If not, I need a standard star reduction, and each .log file needs a corresponding .log reduction that uses the same parameters
    to ensure an accurate match. These should be specified in comp_fnames. If none are supplied, try searching for a


    Arguments:
    ----------
    oname: str
        Template for written files. Applied to

    coords: str
        RA, Dec of target. Must in in a format astropy can understand.

    obsname: str
        Observatory location. See astropy for a list of valid names

    T0: float
        Ephemeris data

    period: float
        Ephemeris data

    SDSS: bool, optional
        If True, I'll do an SDSS lookup for the comparison star magnitudes. If False, use a standard star to calibrate

    std_fname: str, optional
        .log file containing a standard star reduction, to calibrate flux of comparisons

    comp_fnames: list, optional
        list of comparison reductions that match the standard star reduction.

    binsize: int, optional
        Binning of the data.

    myLoc: str, optional
        Working directory. If not supplied, default to current working directory

    ext: list, optional
        Extinction coeffiecients, in order of CCD, mags/airmass

    fnames: list, optional
        List of target reduction files. If not supplied, searches for log files

    std_coords: str, optional
        Ra, Dec of standard star, as a string that Astropy can read

    std_mags: list of float, optional
        Apparent magnitude of the standard in each CCD

    Returns:
    --------
    written_files: list
        List of created .calib files.
    '''

    ## First, find the logfiles we want to use
    if fnames==None:
        fnames = []
        try:
            for file in os.listdir('/'.join([myLoc, 'Reduced_Data'])):
                if file.endswith('.log'):
                    fnames.append('/'.join([myLoc,'Reduced_Data', file]))
        except:
            for file in os.listdir('/'.join([myLoc])):
                if file.endswith('.log'):
                    fnames.append('/'.join([myLoc, file]))
        if len(fnames) == 0:
            printer("  I couldn't find any log files! Stopping...")
            exit()

    if SDSS:
        comp_fnames = [x.replace('.log', '.coords') for x in fnames]
    else:
        if comp_fnames == None:
            "I didn't get any comparison reductions. Please supply these!"
            raise NameError

        # Check we have the same number of comparison reductions as we do target reductions
        if len(comp_fnames) != len(fnames):
            printer("Error! I got a different number of comparison reductions and target reductions!")
            printer("Comparisons:")
            for f in comp_fnames:
                printer(f)
            printer("\nTargets:")
            for f in fnames:
                printer(f)
            raise NameError

    # Writing out
    try:
        os.mkdir('/'.join([myLoc, 'Reduced_Data', 'lightcurves']))
    except: pass

    oname = oname.split('/')
    if oname[0] != 'Reduced_Data':
        oname = ['Reduced_Data'] + oname
    oname = '/'.join(oname)


    # Report the things we're working with
    printer("  Using these log files: ")
    for i, fname in enumerate(fnames):
        printer("    {:2d} - {}".format(i, fname))
    printer('  ')
    printer("  Binning folded data by {}".format(binsize))
    printer("  I'll write out to {}*\n".format(oname))

    #Correct to BMJD
    printer("  Correcting observations from MJD to BMJD (observed from '{}')".format(obsname))

    printer("  Phase folding data for a T0: {:}, period: {:}".format(T0, period))

    # Where are we?
    observatory = coord.EarthLocation.of_site(obsname)
    star_loc = coord.SkyCoord(
        coords,
        unit=(u.hourangle, u.deg), frame='icrs'
    )

    band = ['', 'r',   'g',     'u'   ]
    c    = ['', 'red', 'green', 'blue']
    master = {}
    written_files = []


    # I want a master pdf file with all the nights' lightcurves plotted
    with PdfPages(oname+'_all-nights.pdf') as pdf:
        for fname, refname in zip(fnames, comp_fnames):

            printer("\n----------------------------------------------------------------\n----------------------------------------------------------------\n")
            printer("Calibrating lightcurves for {}".format(fname))
            printer("\n----------------------------------------------------------------\n----------------------------------------------------------------\n")
            data = hcam.hlog.Hlog.from_ascii(fname)

            # Get the apertures of this data set
            aps = data.apnames
            CCDs = [str(i) for i in aps]
            CCDs = sorted(CCDs)

            # If we're in the SDSS field, grab the reference stars' magnitudes from their coords.
            if SDSS:
                printer("  Looking up SDSS magnitudes from the database")
                reference_stars = construct_reference(refname)
            else:
                printer("  Getting comparison star apparent magnitudes, from standard observation")
                reference_stars = get_comparison_magnitudes(std_fname, refname, std_coords=std_coords,
                    comp_coords=coords, std_mags=std_mags, obsname=obsname, ext=ext)

            # Plotting area
            fig, ax = plt.subplots(3, figsize=[12, 8])

            # Loop through the CCDs.
            ### For each CCD, grab the target lightcurve, and the comparisons
            for CCD in CCDs:
                CCD_int = int(CCD)

                printer("  CCD {}".format(CCD))

                # Get this frame's apertures
                ap = aps[CCD]
                # Check that there is more than one aperture -- i.e. if a reference star exists
                if len(ap) == 1:
                    printer("I can't do relative photometry with only one aperture!")
                    printer("!!! Bad log file, '{}' !!!".format(fname))
                    exit()

                # Grab the target data
                target = data.tseries(CCD, '1')

                # First reference star
                reference = data.tseries(CCD, '2')

                # Add up the reference star fluxes
                for comp in ap[2:]:
                    reference = reference + data.tseries(CCD, comp)
                # Take the mean
                reference = reference / len(ap[1:])
                printer("  Instrumental mean counts per frame ({} frames) of {} reference stars: {:.1f}".format(len(reference.y), len(ap[1:]), np.mean(reference.y)))

                ### <reference> is now a mean COUNT of the reference stars for each exposure ###
                ## Calculate their actual mean flux from their apparent magnitudes

                # mags is a list of the relevant reference star magnitudes
                mags = reference_stars[CCD]
                if not SDSS:
                    mags = mags[1:]

                if len(mags) != len(ap[1:]):
                    printer("!!!!!---- len(mags): {} --- len(reference): {}".format(len(mags), len(ap[1:])))

                fluxs = sdss_mag2flux(mags)
                meanFlux = np.mean(fluxs) # Actual FLUX of reference

                printer("  Comparison star magnitudes:".format())
                for m, mag in enumerate(mags):
                    printer("    Star {} -> {:.3f} mag".format(m, mag))
                printer("")
                printer("  Apparent fluxes of the comparison stars:")
                for i, flux in enumerate(fluxs):
                    printer("    Star {} -> {:.3f} mJy".format(i, flux))
                printer('  Mean Flux: {:.3f} mJy\n'.format(meanFlux))
                printer("  mJy per count: {:.3e}".format(meanFlux/np.mean(reference.y)))

                # Conversion of target lightcurve
                reference = reference / meanFlux # Counts/mJy
                ratio = target / reference # mJy

                ratio = tcorrect(ratio, star_loc, obsname)

                if CCD == '1':
                    mintime = np.min(ratio.t)
                    E = calc_E(mintime, T0, period)
                    E = np.rint(E)

                    # The above can be off, if the eclipse isnt the minimum. in/decriment until it's within bounds
                    while T0 + E*period < ratio.t[0]:
                        printer("  !!! Eclipse time not within these data! Attempting to fix...")
                        E += 1
                    while T0 + E*period > ratio.t[-1]:
                        printer("  !!! Eclipse time not within these data! Attempting to fix...")
                        E -= 1

                    printer("  I think that the eclipse spanning from {:.3f} to {:.3f} is cycle number {}".format(
                        ratio.t[0], ratio.t[-1], E)
                    )

                    eclTime = T0 + E*period
                    printer("  The eclipse is then at time {:.3f}".format(eclTime))

                    printer("")


                slice_time = (ratio.t - eclTime) / period
                slice_args = (slice_time < 0.5)  *  (slice_time > -0.5)

                ratio = hcam.hlog.Tseries(
                    slice_time[slice_args],
                    ratio.y[slice_args],
                    ratio.ye[slice_args],
                    ratio.mask[slice_args]
                    )

                printer("  I sliced out {} data from the lightcurve".format(len(ratio.t)))

                # Plotting
                ratio.mplot(ax[CCD_int-1], colour=c[CCD_int])
                ax[CCD_int-1].set_ylabel('Flux, mJy')
                filename = oname
                filename = filename.replace('Reduced_Data', 'Reduced_Data/lightcurves')
                filename = "{}_{}_{}.calib".format(filename, fname.split('/')[-1][:-4], band[CCD_int])
                written_files.append(filename)

                with open(filename, 'w') as f:
                    f.write("# Phase, Flux, Err_Flux\n")
                    for t, y, ye, mask in zip(ratio.t, ratio.y, ratio.ye, ratio.mask):
                        if not mask:
                            f.write("{} {} {}\n".format(t, y, ye))

                # Store the lightcurve in the master dict
                try:
                    # If we find <master> has an entry for this CCD, append it
                    master[CCD] = master[CCD].append(ratio)
                except KeyError:
                    # Otherwise, create a new entry
                    master[CCD] = ratio

                printer("  Finished CCD {}\n".format(CCD))
                del reference
                del target
                del ratio

            ax[0].set_title(fname.split('/')[-1])
            ax[2].set_xlabel('Phase, days')

            plt.tight_layout()
            pdf.savefig()
            plt.close()
    printer("  ")
    plt.close('all')

    return written_files