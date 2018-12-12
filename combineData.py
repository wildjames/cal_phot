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


# ------------------- inputted variables-------------------------------------- #

# ## How do we name the files we write out?
# oname = "MASTER-J2205"

# # Binning
# binsize = 3

# # Star coordinates, observing location [RA, Dec], extinction coefficient
# coords = ['22 05 59.48', '-34 14 33.9']
# obsname = 'lasilla'
# ext = 0.161 ## Mean value over 20 yrs

# # Ephemeris data
# T0 = 58404.1329466701
# period = 0.061286974

# ## List of magnitude corrections ##
# ref_kappa = [np.nan, -25.9314, -26.1793, -24.1461]
# # (First entry is left blank, as CCDs are indexed from 1)

# -------------------\inputted variables-------------------------------------- #

def combineData(oname, coords, obsname, T0, period, ref_kappa=None, SDSS=True, std_fname=None, comp_fnames=None, binsize=10, myLoc='.', ext=0.161, fnames=None, std_coords=None, std_mags=None):
    '''
    oname      - Filename template for writing lightcurve plot and data. Appended with binning factor.
    coords     - RA and DEC of target star. As a string in a format that astropy can interpret.
    obsname    - Observatory name
    T0         - Ephemeris zero point
    period     - Ephemeris period
    ref_kappa  - kappa corrections
    SDSS       - Are we in the SDSS field? If we are, I'll do a lookup for reference magnitudes using coordinates from a '.coords' 
                file of the same name as the logfile (i.e. logfile.log -> logfile.coords), and get the data from there.
    binsize    - Lightcurve binning
    myLoc      - Logfile search directory. Looks here for a directory called 'Reduced_Data', and pulls logfiles from there.
    ext        - Atmospheric extinction, mags/airmass

    Overall goal: Construct a folded lightcurve of an object. 

    Looks in <myLoc> for the directory 'Reduced_Data', and pulls all the logfiles from there for analysis.
    Then, for each observing run, 
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
            print("  I couldn't find any log files! Stopping...")
            exit()
    if comp_fnames == None:
        comp_fnames = fnames
    
    # Check we have the same number of comparison reductions as we do target reductions
    if len(comp_fnames) != len(fnames):
        print("Error! I got a different number of comparison reductions and target reductions!")
        print("Comparisons:")
        for f in comp_fnames:
            print(f)
        print("\nTargets:")
        for f in fnames:
            print(f)
        exit()

    print(comp_fnames)
    
    # Writing out
    try:
        os.mkdir('/'.join([myLoc, 'Reduced_Data', 'lightcurves']))
    except: pass
    oname = oname.split('/')
    if oname[0] != 'Reduced_Data':
        oname = ['Reduced_Data'] + oname
    oname = '/'.join(oname)


    # Report the things we're working with
    print("  Using these log files: ")
    for i, fname in enumerate(fnames):
        print("    {:2d} - {}".format(i, fname))
    print('  ')
    print("  Binning folded data by {}".format(binsize))
    print("  I'll write out to {}*\n".format(oname))

    #Correct to BMJD
    print("  Correcting observations from MJD to BMJD (observed from '{}')".format(obsname))

    print("  Phase folding data for a T0: {:}, period: {:}".format(T0, period))

    # Ephemeris info
    eclFile = 'eclipse_times.txt'
    eclFile = '/'.join([myLoc, eclFile])
    source_key, eclipse_data = read_ecl_file(eclFile)
    
    # Where are we?
    observatory = coord.EarthLocation.of_site(obsname)
    star_loc = coord.SkyCoord(
        coords,
        unit=(u.hourangle, u.deg), frame='icrs'
    )


    if T0 == None:
        T0 = 0.0

    band = ['', 'r',   'g',     'u'   ]
    c    = ['', 'red', 'green', 'blue']
    master = {}
    written_files = []


    # I want a master pdf file with all the nights' lightcurves plotted
    with PdfPages(oname+'_all-nights.pdf') as pdf:
        for fname, refname in zip(fnames, comp_fnames):
            data = hcam.hlog.Hlog.from_ascii(fname)
            
            # Get the apertures of this data set
            aps = data.apnames
            CCDs = [str(i) for i in aps]
            CCDs = sorted(CCDs)

            # If we're in the SDSS field, grab the reference stars' magnitudes from their coords.
            if SDSS:
                refname         = refname.replace('.log', '.coords')
                reference_stars = construct_reference(refname)
            else:
                reference_stars = get_comparison_magnitudes(std_fname, refname, std_coords=std_coords, 
                    comp_coords=coords, std_mags=std_mags, obsname=obsname)

            # Plotting area
            fig, ax = plt.subplots(3, figsize=[12, 8])

            # I want altitude converted to zenith angle. Airmass is roughly constant over 
            # a single eclipse so only do it once to save time.
            obs_T = float(data['1'][0][1])
            obs_T = time.Time(obs_T, format='mjd')
            star_loc_AltAz = star_loc.transform_to(AltAz(obstime=obs_T, location=observatory))
            zenith_angle = 90. - star_loc_AltAz.alt.value
            airmass = 1. / np.cos(np.radians(zenith_angle))
            print("  For the night of {} (observing at {}), calculated altitude of {:.3f}, and airmass of {:.3f}".format(
                fname.split('/')[-1], obs_T.iso, star_loc_AltAz.alt.value, airmass)
            )

            # Loop through the CCDs.
            ### For each CCD, grab the target lightcurve, and the 
            for CCD in CCDs:
                CCD_int = int(CCD)

                # Get this frame's apertures
                ap = aps[CCD]
                # Check that there is more than one aperture -- i.e. if a reference star exists
                if len(ap) == 1:
                    print("  I can't do relative photometry with only one aperture!")
                    print("!!! Bad log file, '{}' !!!".format(fname))
                    exit()

                # Grab the target data
                target = data.tseries(CCD, '1')

                try:
                    # First reference star
                    reference = data.tseries(CCD, '2')
                except:
                    print("Whoops! {} only contains one aperture! I can't do relative photometry with only one!".format(
                        fname
                    ))
                    exit()

                # If we have more than one reference, handle that
                if len(ap) > 2:
                    # Add up the reference star fluxes
                    for comp in ap[2:]:
                        reference = reference + data.tseries(CCD, comp)
                    # Take the mean
                    reference = reference / len(ap[1:])
                
                ### <reference> is now a mean COUNT of the reference stars for each exposure ### 
                ## Calculate their actual mean flux from their apparent magnitudes
                
                # refs is a List of the relevant reference star magnitudes
                mags = reference_stars[CCD]

                fluxs = sdss_mag2flux(mags)
                meanFlux = np.mean(fluxs) # Actual FLUX of reference
                

                print("  CCD {} comparison star magnitudes:".format(CCD))
                for m, mag in enumerate(mags):
                    print("    Comparison {}: {:.3f} mag".format(m, mag))
                print("  ")

                # Conversion of target lightcurve
                reference = reference / meanFlux # Counts/mJy
                ratio = target / reference # mJy

                ratio = tcorrect(ratio, star_loc, obsname)

                if False:
                    # Fold about the period
                    # ratio = ratio.fold(period, t0=T0) ## BUGGED! and not a LTT error?
                    fold_time = (((ratio.t - T0) / period) %1)
                    # fold time domain from -.5 to .5
                    fold_time[fold_time > 0.5] -= 1
                    sorted_args = np.argsort(fold_time)
                    ratio = hcam.hlog.Tseries(
                        fold_time[sorted_args],
                        ratio.y[sorted_args],
                        ratio.ye[sorted_args],
                        ratio.mask[sorted_args]
                        )
                else:
                    ### Cut out an eclipse, without folding it
                    # get E from eclipse data, by interpolating. Only do this for the first CCD, since they should all have the same
                    # timestamps
                    if CCD == '1':
                        E = min([t[1] for t in eclipse_data], key=lambda x:abs(x-np.mean(ratio.t)))
                        E = [t[1] for t in eclipse_data].index(E)
                        E = eclipse_data[E][0]
                        
                        print("  I think that the eclipse spanning from {:.3f} to {:.3f} is cycle number {}".format(
                            ratio.t[0], ratio.t[-1], E)
                        )

                    eclTime = T0 + E*period

                    slice_time = (ratio.t - eclTime) / period
                    slice_args = (slice_time < 0.5)  *  (slice_time > -0.5)

                    ratio = hcam.hlog.Tseries(
                        ratio.t[slice_args] - eclTime,
                        ratio.y[slice_args],
                        ratio.ye[slice_args],
                        ratio.mask[slice_args]
                        )

                # print("  Binning the night of {} by {}".format(fname.split('/')[-1][:-4], binsize))
                # ratio = ratio.bin(binsize)
                
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
                        f.write("{} {} {}\n".format(t, y, ye))

                # Store the lightcurve in the master dict
                try:
                    # If we find <master> has an entry for this CCD, append it
                    master[CCD] = master[CCD].append(ratio)
                except KeyError:
                    # Otherwise, create a new entry
                    master[CCD] = ratio
                

            ax[0].set_title(fname.split('/')[-1])
            ax[2].set_xlabel('Phase, days')

            plt.tight_layout()
            pdf.savefig()
    print("  ")
    plt.close(fig)



    ## <master> is now a dict, containing 3 Tseries objects of the lightcurves.

    # Sort each CCD so all the observations are lined up
    print("  Creating a master lightcurve, by binning the folded, sum lightcurve")
    print("   by (number of eclipses +1): {}".format(len(fnames)+1))
    for CCD in ['1', '2', '3']:
        lightcurve = master[CCD]

        times = lightcurve.t
        sorted_args = np.argsort(times)
        master[CCD] = hcam.hlog.Tseries(
            times[sorted_args],
            lightcurve.y[sorted_args],
            lightcurve.ye[sorted_args],
            lightcurve.mask[sorted_args]
        )

        # Bin the lightcurve by the number of nights
        master[CCD] = master[CCD].bin(len(fnames)+1)

    print("  Done!\n")

    print("  Plotting...")
    if binsize > 1:
        oname += '_bin{:02d}'.format(binsize)
    else:
        oname += "_unbinned"
    # Plotting
    fig, ax = plt.subplots(3, figsize=[24, 8])

    master['1'].mplot(ax[0], colour='red')
    master['2'].mplot(ax[1], colour='green')
    master['3'].mplot(ax[2], colour='blue')

    for x in ax:
        x.set_ylabel('Flux, mJy')
    ax[2].set_xlabel('Phase, days')

    plt.tight_layout()
    plt.savefig(oname+'.pdf')
    # plt.show()
    plt.close('all')
    print('  Saved to {}'.format(oname+'.pdf'))
    print("  Done!\n")

    print("  Writing to files...")

    for i, col in zip(['1', '2', '3'], c[1:]):
        filename = oname+'_{}.calib'.format(col)
        with open(filename, 'w') as f:
            f.write("# Phase, Flux, Err_Flux, Mask\n")
            for t, y, ye, mask in zip(master[i].t, master[i].y, master[i].ye, master[i].mask):
                f.write("{}, {}, {}, {}\n".format(t, y, ye, mask))
        print("  Wrote out {}!".format(filename))

    print("\n  Done!")

    return written_files