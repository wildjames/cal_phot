import copy
import os

import pdb

import hipercam as hcam
import matplotlib.pyplot as plt
import numpy as np
from astropy import coordinates as coord
from astropy import time
from astropy import units as u
from astropy.coordinates import AltAz
from astropy.stats import sigma_clipped_stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit

from .constructReference import construct_reference
from .getEclipseTimes import read_ecl_file, tcorrect
from .logger import printer

## I apologise for this atrocity. This is a project that grew out of some convinience scripts, and by the time it was big enough to be useful to anyone but me, it was a horrid mess. I am very lazy and don't want to clean it up.


# hipercam.FLAGS gives a key for this
FLAGS_TO_IGNORE = [4]


def straight_line(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B

def sdss_mag2flux(mag):
    '''Takes an SDSS magnitude, returns the corresponding flux in [mJy]'''
    alpha = 3631e3

    flux = 10**(-mag/2.5)
    flux*= alpha

    return flux

def sdss_flux2mag(flx):
    '''Takes an flux in mJy and converts it to an SDSS magnitude'''
    alpha = 3631e3

    m = -2.5 * np.log10(flx)
    m += 2.5 * np.log10(alpha)

    return m

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


def extract_data(oname, coords, obsname, T0, period, inst, SDSS,
                 comp_mags=None, myLoc='.', fnames=None,
                 lower_phase=-0.5, upper_phase=0.5,
                 no_calibration=False):
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

    myLoc: str, optional
        Working directory. If not supplied, default to current working directory

    fnames: list, optional
        List of target reduction files. If not supplied, searches for log files

    Returns:
    --------
    written_files: list
        List of created .calib files.
    '''
    printer("\n\n# # # # # # # # # # # # # # # # # # # BEGIN BATCH CALIBRATION # # # # # # # # # # # # # # # # # # #\n\n")

    # Writing out
    lc_dir = os.path.join(myLoc, 'MCMC_LIGHTCURVES')
    try:
        os.mkdir(lc_dir)
    except: pass
    printer("Lightcurves will go in: {}".format(lc_dir))

    figs_dir = os.path.join(myLoc, 'MCMC_LIGHTCURVES', "FIGS")
    try:
        os.mkdir(figs_dir)
    except: pass
    printer("Figures will go in: {}".format(figs_dir))


    # Report the things we're working with
    printer("  Using these log files: ")
    for i, fname in enumerate(fnames):
        printer("    {:2d} - {}".format(i, fname))
    printer('  ')
    printer("  I'll write out to {}*\n".format(oname))

    #Correct to BMJD
    printer("  Correcting observations from MJD to BMJD (observed from '{}')".format(obsname))
    printer("  Phase folding data for a T0: {:}, period: {:}".format(T0, period))

    # Data masking stuff
    FLAG = np.uint32(0)
    for f in FLAGS_TO_IGNORE:
        FLAG = FLAG | f
    if FLAG:
        printer("  Ignoring bad data flags: {}".format(FLAGS_TO_IGNORE))
        printer("List of keys:")
        printer(hcam.FLAGS)

    # Where are we?
    try:
        observatory = coord.EarthLocation.of_site(obsname)
    except:
        lat, lon = obsname.split(',')
        printer("  Attempting to get the earth observatory from latitude and longitude")
        observatory = coord.EarthLocation.from_geodetic(lat=lat, lon=lon)

    star_loc = coord.SkyCoord(
        coords,
        unit=(u.hourangle, u.deg), frame='icrs'
    )

    # I want to know what instrument I'm using, since each has a different number of cameras
    inst = inst.lower()
    if inst == 'uspec':
        nCCD = 1
        bands = ['???']
        c    = ['black']
    elif inst == 'ucam':
        nCCD = 3
        bands = ['r', 'g', 'u']
        c = ['red', 'green', 'blue']
    elif inst == 'hcam':
        nCCD = 5
        bands = ['u', 'g', 'r', 'i', 'z']
        c = ['blue', 'green', 'red', 'magenta', 'black']

    printer("  I'm using the instrument {}, which has {} CCDS in the following order:".format(inst, nCCD))
    for n, col in zip(range(nCCD), c):
        printer("  -> CCD {}: plotted in {}".format(n+1, col))

    written_files = []


    #  Plotting #
    ADU_lightcurves = {fname: [] for fname in fnames}

    print("Making plotting area...", end='')
    plt.ion()
    fig, ax = plt.subplots(nCCD, figsize=[11.69, 8.27], sharex=True)
    # If we only have one CCD, axes still need to be a lists
    if nCCD == 1:
        ax = [ax]

    twinAx = []
    for i, a in enumerate(ax):
        a.set_ylabel('Flux, mJy')

        twinAx.append(a.twinx())
        twinAx[i].set_ylabel('Count Ratio')
        twinAx[i].yaxis.tick_right()

    ax[-1].set_xlabel('Phase, days')
    ax[0].set_title('Waiting for data...')
    fig.tight_layout()

    compFig, compAx = plt.subplots(nCCD, figsize=[11.69,8.27], sharex=True)

    # If we only have one CCD, axes still need to be a lists
    if nCCD == 1:
        compAx = [compAx]

    compFig.tight_layout()
    plt.show()
    print(" Done!")



    # I want a master pdf file with all the nights' lightcurves plotted
    pdfname = os.path.join(figs_dir, oname+"_all_nights.pdf")
    with PdfPages(pdfname) as pdf:
        for fname in fnames:
            printer("\n----------------------------------------------------------------\n----------------------------------------------------------------\n")
            printer("Calibrating lightcurves for {}".format(fname))
            printer("\n----------------------------------------------------------------\n----------------------------------------------------------------\n")

            print("CWD:  {}".format(os.getcwd()))
            print("File: {}".format(fname))
            data = hcam.hlog.Hlog.read(fname)
            if data == {}:
                data = hcam.hlog.Hlog.rulog(fname)
            if data == {}:
                data = hcam.hlog.Hlog.rfits(fname)
            if data == {}:
                raise Exception("Could not properly read in log file, {}".format(fname))

            printer("  Read the data file!")

            # Get the apertures of this data set
            aps = data.apnames
            CCDs = [str(i) for i in aps]
            CCDs = sorted(CCDs)
            if CCDs == []:
                printer("ERROR! No data in the file!")
            printer("  The observations have the following CCDs: {}".format([int(ccd) for ccd in CCDs]))

            printer("  Am I flux calibrating the data? {}".format(not no_calibration))
            if no_calibration:
                printer("\n!!! Not doing flux calibration! Setting reference magnitudes to correspond to a flux=1\n\n")

                # Reference stars is a dict, keyed by the CCD string
                reference_stars = {}
                for CCD in CCDs:
                    mags = []

                    # The comparison star fluxes get added together, so each should have an even share of 1 mJy.
                    individual_flx = 1.0 / len(aps[CCD][1:])
                    for ap in aps[CCD][1:]:
                        mags.append(sdss_flux2mag(individual_flx))
                    reference_stars[CCD] = np.array(mags)

                printer("'Unit' Reference stars have a magnitude of {:.2f}\n".format(sdss_flux2mag(1.0)))

            elif SDSS:
                printer("  Looking up SDSS magnitudes from the database")
                comparison_coord_files = comp_mags[fname]
                reference_stars = construct_reference(comparison_coord_files)

            else:
                reference_stars = {}
                comparisons = comp_mags[fname]
                bands = list(comparisons.keys())

                printer("  For each of these CCDs, I've been given comparison stars of the following magnitudes:")
                for i, (b, comps) in enumerate(comparisons.items()):
                    # I need to capture 'none' strings here, and store them as np.nans.
                    # Later, when I construct the comparison star, these apertures must be ignored!

                    comparison_list = []
                    for comp in comps:
                        try:
                            # Can I float? (we all float down here, georgie...)
                            comp = float(comp)
                            comparison_list.append(comp)
                        except:
                            # If I can't, Ignore me.
                            comparison_list.append(np.nan)
                    # Who doesnt love vectorised calculations?
                    reference_stars[str(i+1)] = np.array(comparison_list)

                printer("  My comparison stars have the following apparent mags:")
                for b, mags in reference_stars.items():
                    printer("    - CCD{}, mags: {}".format(b, mags))
            printer("\n\n")

            for a in ax:
                    a.clear()
                    a.set_ylabel('Flux, mJy')
            ax[-1].set_xlabel('Phase, days')
            ax[0].set_title('Waiting for data...')

            # Loop through the CCDs.
            # # For each CCD, grab the target lightcurve, and the comparisons
            for CCD in CCDs:
                lightcurve_metadata = '# This is data from the file {} CCD {}\n'.format(fname, CCD)

                CCD_int = int(CCD) - 1
                printer("-> CCD {}".format(CCD))

                # Plot the comparison we construct
                compAx[CCD_int].clear()
                compAx[CCD_int].set_title("CCD {}, comparison star".format(CCD))
                compAx[CCD_int].set_ylabel("Counts per frame")

                # Get this frame's apertures
                ap = aps[CCD]
                printer("  This CCD has the apertures: {}".format(ap))
                # Check that there is more than one aperture -- i.e. if a comparison star exists
                if len(ap) <= 1:
                    printer("I can't do relative photometry with only one aperture!")
                    printer("!!! Bad log file, '{}' !!!".format(fname))
                    raise LookupError("Not enough apertures in the log file!", fname)


                # Check for nans in the log files.
                for a in ap:
                    star = data.tseries(CCD, a)
                    if np.any(np.isnan(star.y)):
                        printer("!!! Log file cannot contain nan values! File {}".format(fname))
                        raise ValueError("Log file cannot contain nan values!", fname)


                # Get some data on the """quality""" of the observations
                metadata = '#\n# Reduction info:\n\n'
                to_proc = data[CCD]

                ap_x = [header for header in
                    to_proc.dtype.fields if "x_" in header]
                ap_y = [header for header in
                    to_proc.dtype.fields if "y_" in header]
                ap_fwhm = [header for header in
                    to_proc.dtype.fields if "fwhm_" in header]

                for x_label, y_label, fwhm_label in zip(ap_x, ap_y, ap_fwhm):
                    x_pix_loc = to_proc[x_label].mean()
                    y_pix_loc = to_proc[y_label].mean()
                    fwhm_pix_loc = to_proc[fwhm_label].mean()

                    aperture_number =x_label.replace("x_", "")

                    metadata += '#   Aperture {}\n'.format(aperture_number)
                    metadata += "#     x location: {:.0f}\n".format(x_pix_loc)
                    metadata += "#     y location: {:.0f}\n".format(y_pix_loc)
                    metadata += "#     fwhm:       {:.2f}\n#\n#\n".format(fwhm_pix_loc)
                lightcurve_metadata += metadata


                # Grab the target data
                target = data.tseries(CCD, '1')


                # mags is a list of the relevant comparison star magnitudes.
                # For non-SDSS fields, this is the clipped mean magnitude of each object.
                mags = reference_stars[CCD]
                fluxs = sdss_mag2flux(mags)
                sumFlux = np.nansum(fluxs)
                sumMag = sdss_flux2mag(sumFlux)

                printer("  Comparison star mags: {}".format(mags))
                if no_calibration:
                    lightcurve_metadata += "# No flux calibration being done!!\n"
                    lightcurve_metadata += "# simulated a dummy comparison magnitude of 1.00 mJy\n"
                else:
                    lightcurve_metadata += "# Comparison star mags: {}\n".format(mags)


                # Add up the reference star fluxes
                N_comparisons = 0
                comparison = "Dummy initialiser ( ͡° ͜ʖ ͡°)"
                count_ratios = []
                for a, mag in zip(ap[1:], mags):
                    if np.isnan(mag):
                        printer("  The reference star in ap {} is being ignored!".format(a))
                    else:
                        N_comparisons += 1
                        try:
                            new_comparison = data.tseries(CCD, a)
                            comparison = comparison + new_comparison
                            printer("  The reference star now includes data from aperture {}".format(a))
                        except:
                            comparison = data.tseries(CCD, a)
                            printer("  The comparison was initialised with aperture {}".format(a))
                
                printer("  The 'comparison star' I've construced from {} apertures now has a mean count/frame of {:.3f}".format(N_comparisons, np.mean(comparison.y)))


                # If we have SDSS stars too bright, get their mags from flux calibrating those that arent
                if SDSS and np.any(np.isinf(mags)):
                    printer("\n\n  I have some comparisons that saturated SDSS! Inferring their magnitudes from fainter stars.")

                    lightcurve_metadata += "# Some comparison stars saturated the SDSS image.\n# Their magnitudes were inferred from fainter stars\n"

                    printer("  Collecting fainter stars...")
                    calibComp = None
                    for mag, a in zip(mags, ap[1:]):
                        if np.inf(mag):
                            printer("    Skipping aperture {}, as it is nan".format(a))
                        else:
                            if calibComp is None:
                                calibComp = data.tseries(CCD, a)
                            else:
                                calibComp += data.tseries(CCD, a)

                    if calibComp is None:
                        raise Exception("All comparison stars saturated SDSS! Pick at least one that doesn't!")
                    calibComp_counts = np.mean(calibComp.y)

                    printer("  My non-saturated SDSS stars have a mean count/frame of {:.3f}".format(calibComp_counts))
                    lightcurve_metadata += "# My non-saturated SDSS stars have a mean count/frame of {:.3f}\n".format(calibComp_counts)

                    printer("  My fluxes are {}".format([f for f in fluxs if not np.inf(f)]))
                    printer("    with a sum flux of {:.3f} mJy".format(sumFlux))
                    printer("     and a sum mag of  {:.3f} mag".format(sumMag))

                    lightcurve_metadata += "# My fluxes are {}\n".format([f for f in fluxs if not np.inf(f)])
                    lightcurve_metadata += "#   with a sum flux of {:.3f} mJy\n".format(sumFlux)
                    lightcurve_metadata += "#    and a sum mag of  {:.3f} mag\n".format(sumMag)

                    for i, (mag, a) in enumerate(zip(mags, ap[1:])):
                        if np.isinf(mag):
                            cnts = data.tseries(CCD, a)
                            meanCnts = np.mean(cnts.y)

                            if np.any(np.isnan(cnts.y)):
                                meanCnts = np.nanmean(cnts.y)

                                printer("The file {} has nan counts! That's weird, and you should fix that.".format(fname))
                                printer("I'll continue ignoring the nan, BUT FIX IT!")

                                lightcurve_metadata += "# The file {} has nan counts! That's wierd, and you should fix that.\n".format(fname)
                                lightcurve_metadata += "# I'll continue ignoring the nan, BUT FIX IT!\n"


                            # Calibrated against known SDSS stars. Observed through the same air column since they're the same frame, so no ext. corr.
                            mag = sumMag - 2.5*np.log10(meanCnts/calibComp_counts)
                            mags[i] = mag

                            printer("    Star {} had no SDSS magnitude. Computed a magnitude of {:.3f} from an e- flux of {}".format(a, mag, meanCnts))
                            lightcurve_metadata += "# Star {} had no SDSS magnitude. Computed a magnitude of {:.3f} from an e- flux of {}\n".format(a, mag, meanCnts)

                    printer("\n")

                # # # # # # # # # # # # # # # # # # # #
                # # Conversion of target lightcurve # #
                # # # # # # # # # # # # # # # # # # # #

                # Get the non-saturated fluxes
                fluxs = sdss_mag2flux(mags)
                comparison_flux = np.nansum(fluxs)

                ratio = target / comparison # counts / counts - ratio between target and comp

                printer("\n\n  Correcting data to BMJD time...")
                ratio = tcorrect(ratio, star_loc, obsname)

                # If we're the first CCD, figure out what eclipse cycle we are
                if CCD == '1':
                    meantime = np.mean(ratio.t)
                    E = calc_E(meantime, T0, period)
                    printer("  The mean time of this eclipse is {:.3f}.".format(meantime))
                    printer("  From ephemeris data, I get an eclipse Number,")
                    printer("    E = ({:.3f} - [T0={:.3f}]) / [P={:.5f}]".format(meantime, T0, period))
                    printer("    E = {:.3f}".format(E))

                    E = np.rint(E)
                    # The above can be off, if the eclipse isnt the minimum. in/decriment until it's within bounds
                    while T0 + E*period < ratio.t[0]:
                        printer("    !!! Eclipse time not within these data! Incrimenting E...")
                        E += 1
                    while T0 + E*period > ratio.t[-1]:
                        printer("    !!! Eclipse time not within these data! Decrimenting E...")
                        E -= 1

                    printer("  I think that the eclipse spanning from {:.3f} to {:.3f} is cycle number {}".format(
                        ratio.t[0], ratio.t[-1], E)
                    )

                    eclTime = T0 + E*period
                    printer("  The eclipse is then at time {:.3f}".format(eclTime))
                    printer("")

                # slice out the data between phase -0.5 and 0.5
                printer("  Slicing out data between phase {} and {}".format(lower_phase, upper_phase))
                slice_time = (ratio.t - eclTime) / period
                slice_args = (slice_time < upper_phase)  *  (slice_time > lower_phase)

                ratio = hcam.hlog.Tseries(
                    slice_time[slice_args],
                    ratio.y[slice_args],
                    ratio.ye[slice_args],
                    ratio.mask[slice_args]
                    )

                # Bad data has error = -1
                mask = np.where(ratio.ye != -1)
                ratio = ratio[mask]

                meanRatio = np.mean(ratio.y)

                printer("  I sliced out {} data from the lightcurve about the eclipse.".format(len(ratio.t)))

                # Save the ratio for later
                ADU_lightcurves[fname].append(copy.deepcopy(ratio))

                # Convert the ratio from ADU to mJy
                ratio = ratio * comparison_flux # Scale back up to actual flux.

                # Filter out flags I don't care about.
                ratio.mask = ratio.mask & (~ FLAG)

                #######################################################################################################
                ############ IGNORE ME I'M BORING AND HARD TO READ. WHY READ ANYTHING HARD? JUST TRUST ME. ############
                #######################################################################################################

                #  Reporting  #

                lightcurve_metadata += "# I calculated an eclipse time of {} BMJD, and phase-folded around that\n".format(eclTime)
                lightcurve_metadata += "# with a T0 of {}, and a period of {}, making this eclipse {}.\n".format(T0, period, E)
                lightcurve_metadata += "# I also sliced out the phase {} -> {}!!\n#\n#\n".format(lower_phase, upper_phase)

                printer("  Comparison star apparent SDSS magnitudes:")
                lightcurve_metadata += "# Comparison star apparent SDSS magnitudes:\n"

                for m, mag in enumerate(mags):
                    printer("    Star {} -> {:.3f} mag".format(m, mag))
                    lightcurve_metadata += "#   Star {} -> {:.3f} mag\n".format(m, mag)
                printer("")
                lightcurve_metadata += "#\n#\n#\n"

                printer("  Apparent fluxes of the comparison stars:")
                lightcurve_metadata += "# Apparent fluxes of the comparison stars:\n"

                for i, flux in enumerate(fluxs):
                    printer("    Star {} -> {:.3f} mJy".format(i, flux))
                    lightcurve_metadata += "#   Star {} -> {:.3f} mJy\n".format(i, flux)
                lightcurve_metadata += "#\n"

                printer('  Sum apparent Flux: {:.3f} mJy\n'.format(comparison_flux))

                printer("  Instrumental counts, summed per mean frame ({} frames) of {} comparison stars: {:.1f}".format(
                    len(comparison.y), N_comparisons, np.mean(comparison.y)
                ))

                printer("  Instrumental counts per mean frame ({} frames) of target: {:.1f}".format(
                    len(target.y), np.mean(target.y)
                ))
                printer("  Mean Target/comparison count ratio: {:.3f}".format(meanRatio))
                printer("  Mean target magnitude: {:.3f}".format(sdss_flux2mag(meanRatio * comparison_flux)))

                lightcurve_metadata += '# Sum apparent Flux: {:.3f} mJy\n#\n#\n'.format(comparison_flux)
                lightcurve_metadata += "# Instrumental summed counts per mean frame ({} frames) of {} comparison stars: {:.1f}\n".format(len(comparison.y), len(ap[1:]), np.mean(comparison.y))
                lightcurve_metadata += "# Instrumental counts per mean frame ({} frames) of target: {:.1f}\n#\n".format(len(target.y), np.mean(target.y))
                lightcurve_metadata += "# Mean Target/comparison count ratio: {:.3f}\n".format(meanRatio)
                lightcurve_metadata += "# Mean target magnitude: {:.3f}\n".format(sdss_flux2mag(meanRatio * comparison_flux))

                # Plotting management
                ax[CCD_int].clear()
                if CCD_int == 0:
                    title = os.path.split(fname)[1]

                    ax[0].set_title(title)
                    compAx[0].set_title("{}\nCCD {}, comparison stars, normalised.".format(title, CCD))

                ax[CCD_int].set_ylabel('Flux, mJy')

                # Plot the ratio
                ratio.mplot(ax[CCD_int], colour=c[CCD_int])

                # Scale the right side labels
                twinAx[CCD_int].set_ylim( ax[CCD_int].get_ylim() / comparison_flux )
                # Draw
                fig.canvas.draw_idle()

                compMin =  9e99
                compMax = -9e99
                if len(ap) == 2:
                    # # Plot the mean count flux on the figure -- only used when single aperture, as not as useful as ratios
                    compAx[CCD_int].errorbar(comparison.t, comparison.y, yerr=comparison.ye,
                        label='Mean', color='black', linestyle='', marker='o', capsize=0)
                    try:
                        mean, _, _ = sigma_clipped_stats(comparison.y, maxiters=2, sigma=3)
                    except:
                        mean, _, _ = sigma_clipped_stats(comparison.y, iters=2, sigma=3)
                    compAx[CCD_int].axhline(mean, linestyle='--', color='black')
                    compMin = np.min(comparison.y)
                    compMax = np.max(comparison.y)
                else:
                    # Plot each combination of comparison star ratios, i.e. for 3 comparisons: 2/3, 2/4, 3/4
                    j = 0
                    for i, a in enumerate(ap[1:-1]):
                        first = data.tseries(CCD, a)
                        for b in ap[i+2:]:
                            printer("  -> Plotting ap {}/{}".format(a, b))
                            toPlot = first / data.tseries(CCD, b)

                            # Filter out flags I don't care about.
                            toPlot.mask = toPlot.mask & (~ FLAG)

                            # Apply the mask to the data
                            if np.any(toPlot.mask):
                                mask = np.where(toPlot.mask == 0)
                                printer("  -> {} masked data!".format(np.sum(toPlot.mask != 0)))
                                printer("\nMasked data:")
                                printer(toPlot.mask)
                                printer("\n\n")
                                printer("Flags:")
                                printer(hcam.FLAGS)

                                if np.all(toPlot.mask != 0):
                                    print("ALL DATA ARE MASKED! Stopping...")
                                    exit()

                                toPlot.t  = toPlot.t[mask]
                                toPlot.y  = toPlot.y[mask]
                                toPlot.ye = toPlot.ye[mask]
                                toPlot.mask = toPlot.mask[mask]

                            toPlot.y = toPlot.y / np.mean(toPlot.y)
                            toPlot.y = toPlot.y + (j / 5)
                            j += 1

                            # Get min and max axis limits
                            if np.max(toPlot.y) > compMax:
                                compMax = np.max(toPlot.y)
                            if np.min(toPlot.y) < compMin:
                                compMin = np.min(toPlot.y)

                            # Fit a straight line to the data. Deviations indicate bad comparisons
                            A,B = curve_fit(straight_line, toPlot.t, toPlot.y)[0]
                            fit_X = np.linspace(toPlot.t[0], toPlot.t[-1], 3)
                            fit_Y = straight_line(fit_X, A, B)

                            # iters is depreciated. Try the new version, if that fails do the old version. yay, flexibility!
                            try:
                                mean, _, _ = sigma_clipped_stats(toPlot.y, maxiters=2, sigma=3)
                            except:
                                mean, _, _ = sigma_clipped_stats(toPlot.y, iters=2, sigma=3)
                            compAx[CCD_int].axhline(mean, linestyle='--', color='black')
                            compAx[CCD_int].plot(fit_X, fit_Y, color='red', linestyle=':')
                            compAx[CCD_int].scatter(toPlot.t, toPlot.y,
                                s=10,
                                label="Aperture {}/{} - grad: {:.2f}".format(a, b, A),
                                alpha=0.6
                            )
                            compFig.canvas.draw_idle()

                # Add in legend artist
                compAx[CCD_int].legend()

                pad = 0.05 * (compMax - compMin)
                compAx[CCD_int].set_ylim([compMin-pad, compMax+pad])
                compAx[CCD_int].set_xlim([comparison.t[0], comparison.t[-1]])

                # File handling stuff
                b = bands[CCD_int]
                while b == '???':
                    b = input("What band are these data?: ")
                    if b == '':
                        print("PLEASE ENTER A BAND NAME for:\n{}\n".format(fname))
                        b = "???"

                date = time.Time(eclTime, format='mjd')
                date = date.strftime("%Y-%m-%d@%Hh%Mm")

                filename = oname
                filename = "{}_{}_{}.calib".format(filename, date, b)

                filename = os.path.join(lc_dir, filename)

                # Saving data
                printer("  These data have {} masked points.".format(np.sum(ratio.mask != 0)))
                if np.sum(ratio.mask != 0):
                    printer("\n\n{}\n\n".format(ratio.mask))
                with open(filename, 'w') as f:
                    f.write(lightcurve_metadata)
                    f.write("# Phase, Flux, Err_Flux\n")
                    for t, y, ye, mask in zip(ratio.t, ratio.y, ratio.ye, ratio.mask):
                        if not mask:
                            f.write("{} {} {}\n".format(t, y, ye))

                written_files.append(filename)
                printer("  Wrote data to {}".format(filename))
                printer("  Finished CCD {}\n".format(CCD))

            ax[-1].set_xlabel('Phase, days')

            x_range = [min(ratio.t), max(ratio.t)]
            ax[0].set_xlim(x_range)

            x_range = [min(comparison.t), max(comparison.t)]
            compAx[0].set_xlim(x_range)

            plt.tight_layout()
            plt.subplots_adjust(hspace=0.0)

            fig.canvas.draw_idle()
            compFig.canvas.draw_idle()

            input("\n  Hit enter for next file\r")
            print()
            pdf.savefig(fig)
            pdf.savefig(compFig)
        plt.close(compFig)

        # Plot the ADU lightcurves for the user.
        for a in ax:
                a.clear()
                a.set_ylabel('Flux, ADU')
        ax[-1].set_xlabel('Phase, days')
        ax[0].set_title('Waiting for data...')

        for fname, lightcurves in ADU_lightcurves.items():
            for i, tseries in enumerate(lightcurves):
                print("{} // CCD {}".format(fname, i))
                flx = tseries.y
                phase = tseries.t

                ax[i].step(phase, flx, label=fname)

        for a in ax:
            a.legend()
        ax[0].set_title("ADU Lightcurves of all files")
        plt.tight_layout()
        fig.canvas.draw_idle()
        ADU_name = os.path.join(figs_dir, oname+'_ADU_lightcurves.pdf')
        fig.savefig(ADU_name)

        input("Hit enter to continue... ")

    printer("  ")
    printer("  Saved the plots to {}".format(pdfname))

    plt.close('all')
    plt.ioff()

    return written_files
