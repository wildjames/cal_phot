import hipercam as hcam
import sys
import os
import numpy as np
from astropy import time, units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord


def getKappa(lf, coords, obsname, mags, ext=0.161):
    '''
lf      - Standard star reduced data in a logfile from the hipercam pipeline
mags    - [<r' magnitude>, <g' band magnitude>, <u' band magnitude>] of standard
coords  - RA and DEC of standard star, in hours and degrees, respectively.
obsname - name of observatory
ext     - optional parameter, if you want to use a different extinction coeff.


Takes a logfile containing the reduced lightcurve of a standard star, and a list of its magnitudes.
Uses the observatory and coordinates of the star to get its airmass at the start of the run, and 
assumes it's rouchly constant over the observations.

From these, calculate the instrumental magnitude of the standard, and get the correction factor, kappa,
in each band. Returns these as a list.
'''

    # Double check our data are the right format
    extinction_coefficient = float(ext)
    mags = [float(i) for i in mags]

    # Just checking
    if not os.path.isfile(lf):
        print("Couldn't find the log file:\n{}".format(lf))
        exit()

    # Where are we?
    observatory = EarthLocation.of_site(obsname)
    star_loc = SkyCoord(
        coords,
        unit=(u.hour, u.deg), frame='icrs'
    )

    # Storage list. 
    mags = [np.nan] + mags

    # grab the data
    data = hcam.hlog.Hlog.from_ascii(lf)

    ### Airmass ###
    # Get the time of the observations.
    # ----> Assume there's a negligible change in airmass over the run.
    T = float(data['1'][0][1])
    T = time.Time(T, format='mjd')
    date = T.to_datetime()
    # Calculate the altitude at that time
    star_loc_AltAz = star_loc.transform_to(AltAz(obstime=T, location=observatory))
    zenith_angle = 90. - star_loc_AltAz.alt.value
    
    # Calculate the airmass
    airmass = 1. / np.cos(np.radians(zenith_angle))
    print("  For the observations on {}, calculated z of {:.3f}, and airmass of {:.3f}.\n".format(
        date.strftime("%Y-%m-%d"), zenith_angle, airmass)
    )

    kappas = [np.nan]
    # Get the fluxes from each CCD
    for CCD in ['1', '2', '3']:
        fluxs = []
        star = data.tseries(CCD, '1')
        for counts, line in zip(star.y, data[CCD]):
            # Third column contains exposure time
            flux = counts / float(line[3])
            fluxs.append(flux)

        flux = np.mean(fluxs)
        # print("For CCD{}, found a mean flux of {:.2f} counts/s".format(CCD, flux))
        
        # Instrumental magnitude
        inst_mag = -2.5*np.log10(flux) - (extinction_coefficient*airmass)

        # Correction factor for this band
        kappa = inst_mag - mags[int(CCD)]

        # print("This corresponds to a correction of {:.3} magnitudes.\n".format(kappa))
        kappas.append(kappa)

    return kappas