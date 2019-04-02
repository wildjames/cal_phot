import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np

import os
import copy

from astropy import time, coordinates as coord, units as u
from astropy.coordinates import AltAz
from astropy.stats import sigma_clipped_stats

import hipercam as hcam
from constructReference import construct_reference, get_comparison_magnitudes
from getEclipseTimes import read_ecl_file, tcorrect
from logger import printer

def sdss_mag2flux(mag):
    '''Takes an SDSS magnitude, returns the corresponding flux in [mJy]'''
    alpha = 3631e3

    flux = 10**(-mag/2.5)
    flux*= alpha

    return flux

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


def combineData(oname, coords, obsname, T0, period, inst='ucam', SDSS=True, std_fname=None, comp_fnames=None,
                myLoc='.', ext=None, fnames=None, std_coords=None, std_mags=None):
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
        os.mkdir('/'.join([myLoc, 'lightcurves']))
    except: pass
    try:
        os.mkdir('/'.join([myLoc, 'figs']))
    except: pass

    oname = oname.split('/')
    if oname[0] != myLoc:
        oname = [myLoc, 'lightcurves'] + oname
    oname = '/'.join(oname)


    # Report the things we're working with
    printer("  Using these log files: ")
    for i, fname in enumerate(fnames):
        printer("    {:2d} - {}".format(i, fname))
    printer('  ')
    printer("  I'll write out to {}*\n".format(oname))

    #Correct to BMJD
    printer("  Correcting observations from MJD to BMJD (observed from '{}')".format(obsname))

    printer("  Phase folding data for a T0: {:}, period: {:}".format(T0, period))

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

    inst = inst.lower()
    if inst == 'uspec':
        nCCD = 1
        band = ['???']
        c    = ['black']
    elif inst == 'ucam':
        nCCD = 3
        band = ['r', 'g', 'u']
        c = ['red', 'green', 'blue']
    elif inst == 'hcam':
        nCCD = 5
        band = ['u', 'g', 'r', 'i', 'z'] #TODO: verify this!!
        c = ['blue', 'green', 'red', 'magenta', 'black']

    printer("  I'm using the instrument {}, which has {} CCDS in the following order:".format(inst, nCCD))
    for n, b, col in zip(range(nCCD), band, c):
        printer("  -> CCD {}: {} band, plotted in {}".format(n+1, b, col))

    master = {}
    written_files = []


    # Plotting area
    print("Making plotting area...", end='')
    plt.ion()
    fig, ax = plt.subplots(nCCD, figsize=[11.69,8.27], sharex=True)

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
        ax = [ax]

    compFig.tight_layout()
    plt.show()
    print(" Done!")

    # I want a master pdf file with all the nights' lightcurves plotted
    pdfname = '/'.join([myLoc, 'figs', 'all_nights.pdf'])
    with PdfPages(pdfname) as pdf:
        for fname, refname in zip(fnames, comp_fnames):
            printer("\n----------------------------------------------------------------\n----------------------------------------------------------------\n")
            printer("Calibrating lightcurves for {}".format(fname))
            printer("\n----------------------------------------------------------------\n----------------------------------------------------------------\n")
            try:
                data = hcam.hlog.Hlog.from_ascii(fname)
                if data == {}:
                    raise Exception
            except:
                data = hcam.hlog.Hlog.from_ulog(fname)

            # Get the apertures of this data set
            aps = data.apnames
            CCDs = [str(i) for i in aps]
            CCDs = sorted(CCDs)
            if CCDs == []:
                printer("ERROR! No data in the file!")
            printer("  The observations have the following CCDs: {}".format([int(ccd) for ccd in CCDs]))

            # If we're in the SDSS field, grab the reference stars' magnitudes from their coords.
            if SDSS:
                printer("  Looking up SDSS magnitudes from the database")
                reference_stars = construct_reference(refname)
            else:
                printer("  Getting comparison star apparent magnitudes, from standard observation")
                reference_stars = get_comparison_magnitudes(std_fname, refname, std_coords=std_coords,
                    comp_coords=coords, std_mags=std_mags, obsname=obsname, ext=ext)

            for a in ax:
                    a.clear()
                    a.set_ylabel('Flux, mJy')
            ax[-1].set_xlabel('Phase, days')
            ax[0].set_title('Waiting for data...')

            # Loop through the CCDs.
            ### For each CCD, grab the target lightcurve, and the comparisons
            for CCD in CCDs:
                CCD_int = int(CCD) - 1

                printer("-> CCD {}".format(CCD))

                # Get this frame's apertures
                ap = aps[CCD]
                printer("  This CCD has the apertures: {}".format(ap))
                # Check that there is more than one aperture -- i.e. if a comparison star exists
                if len(ap) == 1:
                    printer("I can't do relative photometry with only one aperture!")
                    printer("!!! Bad log file, '{}' !!!".format(fname))
                    exit()

                # Grab the target data
                target = data.tseries(CCD, '1')

                # First comparison star
                comparison = data.tseries(CCD, '2')
                printer("  Got the reference star in aperture 2")

                # Plot the comparison we construct
                compAx[CCD_int].clear()
                compAx[CCD_int].set_title("CCD {}, comparison star".format(CCD))
                compAx[CCD_int].set_ylabel("Counts per frame")

                # Add up the reference star fluxes
                N = 1
                for a in ap[2:]:
                    N += 1
                    comparison = comparison + data.tseries(CCD, a)
                    printer("  The reference star now includes data from aperture {}".format(a))

                # Take the mean
                comparison = comparison / N

                ### <comparison> is now a mean COUNT of the comparison stars for each exposure ###
                ## Calculate their actual mean flux from their apparent magnitudes

                # mags is a list of the relevant comparison star magnitudes.
                # For non-SDSS fields, this is the clipped mean magnitude of each object.
                mags = reference_stars[CCD]


                if len(mags) != N:
                    printer("  Target reduction filename:    {}".format(fname))
                    printer("  Standard reduction filename:  {}".format(refname))
                    printer("!!!!! Number of comparison magnitudes in standard star reduction: {}".format(len(mags)))
                    printer("!!!!! Number of comparison stars in target reduction: {}".format(len(ap[1:])))
                    input("Hit <Enter> to continue")

                fluxs = sdss_mag2flux(mags)
                # Don't take the clipped mean here, as the meanFlux is the mean, mean flux of our comparison stars,
                #
                meanFlux  = np.mean(fluxs) # Actual FLUX of constructed comparison star

                ## Reporting
                printer("  Comparison star apparent SDSS magnitudes:".format())
                for m, mag in enumerate(mags):
                    printer("    Star {} -> {:.3f} mag".format(m, mag))
                printer("")
                printer("  Apparent fluxes of the comparison stars:")
                for i, flux in enumerate(fluxs):
                    printer("    Star {} -> {:.3f} mJy".format(i, flux))
                printer('  Mean apparent Flux: {:.3f} mJy\n'.format(meanFlux))
                printer("  Instrumental mean counts per frame ({} frames) of {} comparison stars: {:.1f}".format(
                    len(comparison.y), len(ap[1:]), np.mean(comparison.y)
                ))
                printer("  mJy per count: {:.3e}".format(meanFlux/np.mean(comparison.y)))


                # Conversion of target lightcurve
                cnt_per_flx = comparison / meanFlux # Counts/mJy
                ratio = target / cnt_per_flx # mJy

                printer("  Correcting data to BMJD time...")
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

                printer("  I sliced out {} data from the lightcurve about the eclipse.".format(len(ratio.t)))

                # Plotting management
                ax[CCD_int].clear()
                if CCD_int == 0:
                    ax[0].set_title(fname.split('/')[-1])
                    compAx[0].set_title("{}\nCCD {}, comparison star".format(fname.split('/')[-1], CCD))
                ax[CCD_int].set_ylabel('Flux, mJy')

                # Scale the right labels
                twinAx[CCD_int].set_ylim( ax[CCD_int].get_ylim() / meanFlux )

                # Plot the ratio
                ratio.mplot(ax[CCD_int], colour=c[CCD_int])

                # # Plot the mean count flux on the figure -- only used when single aperture, as not as useful as ratios
                if len(ap) == 2:
                    compAx[CCD_int].scatter(comparison.t, comparison.y, s=5, label='Mean', color='black')
                    compAx[CCD_int].axhline(np.mean(comparison.y), linestyle='--', color='black')
                else:
                    # Plot each combination of comparison star ratios, i.e. for 3 comparisons: 2/3, 2/4, 3/4
                    for i, a in enumerate(ap[1:-1]):
                        first = data.tseries(CCD, a)
                        for b in ap[i+2:]:
                            print("Plotting ap {}/{}".format(a, b))
                            toPlot = first / data.tseries(CCD, b)

                            mean, _, _ = sigma_clipped_stats(toPlot.y, maxiters=2, sigma=3)
                            compAx[CCD_int].axhline(mean, linestyle='--', color='black')
                            compAx[CCD_int].scatter(toPlot.t, toPlot.y,
                                s=10,
                                label="Aperture {}/{}".format(a, b),
                                alpha=0.6,
                                color='red'
                            )

                # Add in legend artist
                compAx[CCD_int].legend()

                # recompute the ax.dataLim
                compAx[CCD_int].relim()
                # update ax.viewLim using the new dataLim
                compAx[CCD_int].autoscale_view()

                # File handling stuff
                b = band[CCD_int]
                if b == '???':
                    b = ''
                else:
                    b = '_'+b

                filename = oname
                filename = "{}_{}{}.calib".format(filename, fname.split('/')[-1][:-4], b)

                # Saving data
                with open(filename, 'w') as f:
                    f.write("# Phase, Flux, Err_Flux\n")
                    for t, y, ye, mask in zip(ratio.t, ratio.y, ratio.ye, ratio.mask):
                        if not mask:
                            f.write("{} {} {}\n".format(t, y, ye))

                written_files.append(filename)

                # Store the lightcurve in the master dict
                try:
                    # If we find <master> has an entry for this CCD, append it
                    master[CCD] = master[CCD].append(ratio)
                    printer("  Appended CCD {} to the data dict.".format(CCD))
                except KeyError:
                    # Otherwise, create a new entry
                    printer("  Adding {} to the data dict...".format(CCD))
                    master[CCD] = ratio


                printer("  Finished CCD {}\n".format(CCD))

            ax[-1].set_xlabel('Phase, days')

            x_range = [min(ratio.t), max(ratio.t)]
            ax[0].set_xlim(x_range)

            x_range = [min(comparison.t), max(comparison.t)]
            compAx[0].set_xlim(x_range)

            plt.tight_layout()

            fig.canvas.draw_idle()
            compFig.canvas.draw_idle()

            input("\n  Hit enter for next file\r")
            print()
            pdf.savefig(fig)
            pdf.savefig(compFig)
    printer("  ")
    plt.close('all')
    plt.ioff()

    return written_files