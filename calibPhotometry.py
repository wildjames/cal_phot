import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import json
import os
from astropy import time, coordinates as coord, units as u

import hipercam as hcam
from constructReference import construct_reference

'''
fname: 
    This is the .log file that is produced when reducing data with the HiPERCAM 
    pipeline. I pull the photometry outputs from here and process them
refFname:
    This is a list of coordinates, in plaintext in the following format:
        <CCD 1 reference 1 RA> <CCD 1 reference 1 Dec>
        <CCD 1 reference 2 RA> <CCD 1 reference 2 Dec>
        
        <CCD 2 reference 1 RA> <CCD 2 reference 1 Dec>
        <CCD 2 reference 2 RA> <CCD 2 reference 2 Dec>
        <CCD 2 reference 3 RA> <CCD 2 reference 3 Dec>
        
        <CCD 3 reference 1 RA> <CCD 3 reference 1 Dec>
        <CCD 3 reference 2 RA> <CCD 3 reference 2 Dec>
    # are comments, if necessary, and the next CCD is denoted by a blank line. If there are 
    more than three blocs, only the first three are read in.
oname:
    This is a template name used to write out the .pdf and the .calib files

Creates a .pdf with the calibrated lightcurve of aperture 1 (from the logfile supplied), using the 
other apertures as reference. Requires calibration against the true magnitude of the reference 
stars, which are pulled from SDSS using the coordinates given in the file <refFname>. Also writes 
'<oname>.calib', which contains all the reduced data in a csv with the headers:
    MJD,   r_Mag,   Err_r_Mag,   g_Mag,   Err_g_Mag,   u_Mag,   Err_u_Mag

'''

# -------------------User inputted variables----------------------------------------- #

# Orbital Period
T0 = 57889.097660085354
period = 0.375538871
# binning
binsize = 10

# Observing location
obsname = input("What observatory were the observations made: ")

## How do we name the files we write out?
oname = input('Please enter a target name\n(will be used for writing out files in the "Reduced_Data" directory): ')

# -------------------\User inputted variables----------------------------------------- #


oname = oname.split('/')
if oname[0] != 'Reduced_Data':
    oname = ['Reduced_Data'] + oname
oname = '/'.join(oname)
print("I'll write out to {}".format(oname))


def sdss_mag2flux(mag):
    '''Takes an SDSS magnitude, and returns the corresponding flux in [mJy]'''
    alpha = 3631.e3

    flux = 10**(-mag/2.5)
    flux*= alpha

    return flux

# Setup useful lists
band = ['', 'r',   'g',     'u'   ]
c    = ['', 'red', 'green', 'blue']


# The meat of it

## First, find the logfiles we want to use
fnames = []
try:
    for file in os.listdir('Reduced_Data/'):
        if file.endswith('.log'):
            fnames.append('/'.join(['.','Reduced_Data', file]))
except:
        print("I had trouble finding the 'Reduced_Data' folder - does it exist?")
if len(fnames) == 0:
        print("I couldn't find any log files! Stopping...")
        exit()
# List the files we found
print("Found these log files: ")
for i, fname in enumerate(fnames):
    print("{:2d} - {}".format(i, fname))
print('')

print("Found the following coordinate lists:")
for i, fname in enumerate(fnames):
    # the cood_list files should be in teh same directory as the log file they refer to
    refFname = fname.replace('.log', '.coords')
    if os.path.isfile(refFname):
        print('{:2d} - {}'.format(i, refFname))
print('')

if input('OK to proceed? y/n: ').lower() != 'y':
    exit()




pdfname = oname+'.pdf'
for fname in fnames:
    data = hcam.hlog.Hlog.from_ascii(fname)
    reference_stars = construct_reference(
        fname.replace('.log', '.coords')
    )

    # Get the apertures for all the CCDs
    aps = data.apnames

    for CCD in aps:
        CCD_int = int(CCD)
        # Get the apertures for this CCD
        ap = aps[CCD]
        # Check that there is more than one aperture -- i.e. if a reference star exists
        if len(ap) == 1:
            print("I can't do relative photometry with only one aperture!")
            continue


        # Grab the target data
        target = data.tseries(CCD, '1')

        # First reference star
        reference = data.tseries(CCD, '2')

        # If we have more than one reference, handle that
        if len(ap) > 1:
            for comp in ap[2:]:
                reference = reference + data.tseries(CCD, comp)
            # Take the mean
            reference = reference / len(ap[1:])
        
        ### <reference> is now a mean COUNT of the reference stars ### 
        ## Calculate their actual mean flux...
        
        # List of the relevant reference star data
        refs = reference_stars[CCD]

        # Construct an array of the reference magnitudes, with some super opaque list comprehension.
        # Just trust me on this one
        mags = np.array(
            [ float(
                refs[comp][ band[CCD_int] ]
            )
                for comp in refs ]
        )

        fluxs = sdss_mag2flux(mags)
        meanFlux = np.mean(fluxs) # Actual FLUX of reference

        # Conversion array
        reference = reference / meanFlux # Counts/mJy
        ratio = target / reference # mJy

        ## Check out master dict exists
        try:
            master # Fails if master doesn't exist
        except:
            # <master> doesn't exist, so initialise it as the reference we just calculated
            master = {}
            master[CCD] = ratio
        else:
            ## Check we can append
            try:
                # If we find <master> has an entry for this CCD, append it
                master[CCD] = master[CCD].append(ratio)
            except:
                # Otherwise, create a new entry
                master[CCD] = ratio

## I need to correct MJD to BMJD. Use astropy for this...
# I need the observation coordinates. The reference stars are close, so I'll use them.
coords = [refs['2']['ra'], refs['2']['dec']]


## <master> is now a dict, containing 3 Tseries objects of the lightcurves.
#Correct to BMJD
print("Correcting observations from MJD to BMJD (observed from '{}'), for the coordinates:".format(obsname))
print("RA: {}\nDec: {}".format(coords[0], coords[1]))

# Where are we?
observatory = coord.EarthLocation.of_site(obsname)
star_loc = coord.SkyCoord(
    ra=coords[0], dec=coords[1],
    unit=(u.hourangle, u.deg), frame='icrs'
)

print("And folding data for a T0: {:.6f}, period: {:.6f}".format(T0, period))
# Apply the correction to each CCD
for CCD in ['1', '2', '3']:
    times = time.Time(master[CCD].t,
        format='mjd', scale='utc', location=observatory
    )
    ltt = times.light_travel_time(star_loc)
    times = times.tdb + ltt
    times = np.array(
        [t.value for t in times]
    )

    # # Manually fold the times
    # times = ((times - T0)/period) % 1.
    # times[times > 0.5] -= 1.

    # sorted_args = np.argsort(times)

    # master[CCD] = hcam.hlog.Tseries(
    #     times[sorted_args],
    #     master[CCD].y[sorted_args],
    #     master[CCD].ye[sorted_args],
    #     master[CCD].mask[sorted_args]
    # )

    master[CCD] = master[CCD].fold(period, t0=T0)
    # master[CCD] = master[CCD].remove_outliers()
    master[CCD] = master[CCD].bin(binsize)

print("Done!\n")

print("Plotting...")
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
ax[-1].set_xlabel('Phase, days')

plt.tight_layout()
plt.savefig(oname+'.pdf')
print('Saved to {}'.format(oname+'.pdf'))
print("Done!\n")

print("Writing to file")
with open(oname+'.calib', 'w') as f:
    for i, col in zip(['1', '2', '3'], c[1:]):
        f.write("# {} band observations\n".format(col))
        f.write("# Phase, Flux, Err_Flux, Mask\n")
        for t, y, ye, mask in zip(master[i].t, master[i].y, master[i].ye, master[i].mask):
            f.write("{}, {}, {}, {}\n".format(t, y, ye, mask))
print("Wrote to {}".format(oname+'.calib'))
print("Done!")
