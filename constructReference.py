import json
import requests
import os
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt

from astropy import time, coordinates as coord, units as u
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import AltAz

import hipercam as hcam
from logger import printer

'''

This script is going to take a list of co-ordinates, and use that to query SDSS for the magnitudes
of those stars.

SDSS REST format:
    GET or POST /radialSearch
    Prameters      - Expected Values
    Ra             - Right Ascention in degrees
    Dec            - Declination in degrees
    radius         - Radius in arcminutes
    format         - output file format. E.g. csv,html,xml
    whichway       - Specify Equitorial or Galactic
    uband          - Specify comma seperated range e.g 0,20. This is an optional parameter for SDSS U band.
    gband          - Specify comma seperated range e.g 0,20. This is an optional parameter for SDSS G band.
    rband          - Specify comma seperated range e.g 0,20. This is an optional parameter for SDSS R band.
    iband          - Specify comma seperated range e.g 0,20. This is an optional parameter for SDSS I band.
    zband          - Specify comma seperated range e.g 0,20. This is an optional parameter for SDSS Z band.
    whichquery     - imaging or spectra

    Example:
    http://skyserver.sdss.org/dr12/SkyserverWS/SearchTools/
    RadialSearch?ra=258.2&dec=64&radius=4.1&whichway=equitorial&
    limit=10&format=json&fp=none&uband=0,17&gband=0,15&whichquery=imaging

It will construct a dict of the magnitudes from the result and return it.

'''

def robust_mag(cps):
    '''Converts a list of fluxes

    '''
    mean, median, sigma = sigma_clipped_stats(cps, iters=2, sigma=3)
    return -2.5*np.log10(mean)


def deg2arcsec(inp, ra):
    sign = '+'
    if inp < 0:
        inp = -inp
        sign = '-'

    if ra:
        hh = (24./360.) * inp
    else:
        hh = inp
    mm = (hh - int(hh)) * 60
    ss = (mm - int(mm)) * 60

    hh = int(hh)
    mm = int(mm)

    output = '{}{:d}:{:d}:{:.2f}'.format(sign, hh, mm, ss)

    return output


def construct_reference(fetchFname):
    '''
    Queries the SDSS database for the magnitudes of the stars contained in <fetchFname>, returns them as a dict of lists.

    mags = {
        '1': [ap1, ap2, ap3],      # r' band
        '2': [ap1, ap2],           # g' band
        '3': [ap1, ap2, ap3, ap4]  # u' band
    }

    <fetchFname> is formatted as follows:
        "CCD 1 reference 1 RA" "CCD 1 reference 1 Dec"
        "CCD 1 reference 2 RA" "CCD 1 reference 2 Dec"

        "CCD 2 reference 1 RA" "CCD 2 reference 1 Dec"
        "CCD 2 reference 2 RA" "CCD 2 reference 2 Dec"
        "CCD 2 reference 3 RA" "CCD 2 reference 3 Dec"

        "CCD 3 reference 1 RA" "CCD 3 reference 1 Dec"
        "CCD 3 reference 2 RA" "CCD 3 reference 2 Dec"


    Arguments:
    ----------
    fetchFname: string
        File containing the RA and Dec of the stars, in the format above


    Returns:
    --------
    mags: dict
        Dictionary containing the magnitudes of the stars in the relevant band.
    '''

    printer("\n\n--- Getting reference SDSS magnitudes from '{:s}' ---".format(fetchFname.split('/')[-1]))

    CCDs = ['1', '2', '3']
    radius = '0.1' # search radius, arcseconds

    fetchme = {
        "1": [],
        "2": [],
        "3": []
    }

    if not os.path.isfile(fetchFname):
        printer("The file I was passed, {}, does not exist!".format(fetchFname))
        # Check that we;ve not been given a directory with a 'coord_list' file in it:
        test = fetchFname.split('/')
        test.append('coord_list.coords')
        test = '/'.join(test)
        if os.path.isfile(test):
            fetchFname = test
            printer("I did, however, find a file called {}. Using that instead...".format(fetchFname))
        # If not, then we have no file. Create a template.
        else:
            with open(fetchFname, 'w') as f:
                f.write('<CCD 1 reference 1 RA> <CCD 1 reference 1 Dec>\n')
                f.write("<CCD 1 reference 2 RA> <CCD 1 reference 2 Dec>\n")
                f.write("\n")
                f.write("<CCD 2 reference 1 RA> <CCD 2 reference 1 Dec>\n")
                f.write("<CCD 2 reference 2 RA> <CCD 2 reference 2 Dec>\n")
                f.write("<CCD 2 reference 3 RA> <CCD 2 reference 3 Dec>\n")
                f.write("\n")
                f.write("<CCD 3 reference 1 RA> <CCD 3 reference 1 Dec>\n")
                f.write("<CCD 3 reference 2 RA> <CCD 3 reference 2 Dec>\n")
            printer("Couldn't find that co-ordinate list! I created a template for you at the location you gave me, {}".format(fetchFname))
            raise FileNotFoundError

    printer("Getting SDSS magnitudes for the coordinates found in {}".format(fetchFname))

    with open(fetchFname) as f:
        x = 1
        for line in f:
            if line[0] == '#':
                # print(line.strip())
                pass
            elif len(line.split()) != 2:
                x += 1
            else:
                fetchme[str(x)].append(line.split())

    toWrite = {
        '1':[],
        '2':[],
        '3':[]
    }

    bands = ['', 'r', 'g', 'u']

    for CCD in CCDs:
        printer('-> CCD {}'.format(CCD))
        # Grab the list of coordinates we want to query
        coords = fetchme[CCD]

        for i, coord in enumerate(coords):
            ra, dec = coord

            printer('    Searching -> RA, Dec: {}, {}'.format(ra, dec))

            # Construct the URL we're gonna post. First define what DB we want to search
            url  = 'http://skyserver.sdss.org/dr14/SkyserverWS/SearchTools/RadialSearch?'
            url += 'ra={}&dec={}&'.format(ra, dec)
            # I'm using a radial search, this is the readius of that search
            url += 'radius={}&'.format(radius)
            # Which coord system are we using. Can't imagine I'll want galactic...
            url += 'whichway=equitorial&'
            # If I'm getting more than 5 results, I've probably picked a crappy reference anyway...
            url += 'limit=5&'
            url += 'format=json&'
            url += 'whichquery=imaging'

            resp = requests.post(url)

            results = resp.json()[0]['Rows']
            if len(results) >= 5:
                printer('You got a lot of results from the SDSS query! Choose from the following VERY carefully.')
            if len(results) > 1:
                printer("--------------------------------------------\nMore than one object found at that location!")
                # Get the user to pick one:
                for m, line in enumerate(results):
                    printer("{}\n  RA: {}, Dec: {}\n  u: {}\n  g: {}\n  r: {}".format(
                        m, line['ra'], line['dec'],
                        line['u'], line['g'], line['r']
                        )
                    )
                n = input("Which object to use?: ")
                printer("Chose object {}".format(n), terminal=False)
                target = results[int(n)]
                printer('--------------------------------------------')
            elif len(results) == 1:
                target = results[0]
                ra = deg2arcsec(target['ra'], ra=True)
                dec = deg2arcsec(target['dec'], ra=False)
                printer("    Found one result:\n      ra: {}, dec: {}\n        u: {}\n        g: {}\n        r: {}".format(
                        ra, dec,
                        target['u'], target['g'], target['r']
                        )
                    )
            else:
                printer('ERROR! Found no targets at the location: RA: {}, Dec: {}'.format(target['ra'], target['dec']))
                printer('Try broadening the search radius in this script (was {}),'.format(radius))
                printer('and make sure that your targets are definitely in the SDSS field!')
                raise LookupError

            # pprint(target)
            toWrite[CCD].append(
                target[ bands[int(CCD)] ]
            )

        toWrite[CCD] = np.array(toWrite[CCD])

    printer("Done!\n")
    return toWrite

def get_instrumental_mags(data, coords=None, obsname=None, ext=None):
    '''
    Takes a hipercam data object, and exctracts the instrumental magnitude of each aperture in each CCD

    If Coords and an observatory are supplied, also correct for airmass, using supplied extinction coeffiecients


    Arguments:
    ----------
    data: hipercam.Hlog
        The data to analyse. tseries will be extracted from here.

    coords: str
        Optional. Ra, Dec of the data. Must be readable by Astropy.

    obsname: str
        Optional. Observing location name of the data.

    ext: list-like
        Optional. Extinction corrections to apply, given in CCD order. i.e. [<CCD1 ext>, <CCD2 ext>, ...]


    Returns:
    --------
    all_mags: dict
        Dict, with the keys corresponding to the CCD numbers. Each entry is a numpy array of the
        instrumental magnitudes, in the order they're found in the aperture list.
    '''
    printer("------- Getting instrumental magnitude -------")

    if coords != None and obsname != None:
        printer("  I'm correcting for airmass, using the following:")
        printer("     Extinction: {} mags/airmass".format(ext))
        printer("        Ra, Dec: {}".format(coords))
        printer("    Observatory: {}".format(obsname))

        # Where are we?
        observatory = coord.EarthLocation.of_site(obsname)
        star_loc = coord.SkyCoord(
            coords,
            unit=(u.hourangle, u.deg), frame='icrs')

        # I want altitude converted to zenith angle. Airmass is roughly constant over
        # a single eclipse so only do it once to save time.
        obs_T = float(data['1'][0][1])
        obs_T = time.Time(obs_T, format='mjd')
        star_loc_AltAz = star_loc.transform_to(AltAz(obstime=obs_T, location=observatory))
        zenith_angle = (np.pi/2)*u.rad - (star_loc_AltAz.alt.rad * u.rad)
        airmass = 1. / np.cos(zenith_angle)
        printer("  For the observations at {}, calculated altitude of {:.3f}, and airmass of {:.3f}\n".format(
            obs_T.iso, star_loc_AltAz.alt.value, airmass))
    else:
        printer("  No coordinates or observatory provided, setting airmass to 0.0")
        airmass = 0.0


    all_mags = {}
    aps = data.apnames
    CCDs = [str(i+1) for i, key in enumerate(aps)]

    if ext == None:
        ext = [0.0 for i in CCDs]

    for CCD in CCDs:
        # Get this frame's apertures
        ap = sorted(aps[CCD])

        ex = ext[int(CCD)-1]

        star = data.tseries(CCD, '1')

        # star counts/s
        fl = star.y / data[CCD]['Exptim']

        # Calculate the mean apparent magnitude of the star above the atmosphere
        mag = robust_mag(fl)
        # star magnitudes
        mags = [mag]

        # mean, median, sigma = sigma_clipped_stats(fl, iters=2, sigma=3)
        # regMean = np.mean(fl)
        # printer("CCD {}, Aperture {} clipped mean count flux: {:.3f}".format(CCD, '1', mean))

        # fig, ax = plt.subplots(figsize=[8,6])
        # ax.scatter([i for i, val in enumerate(fl)], fl, color='black')
        # # Poissonian error only
        # ax.errorbar([i for i, val in enumerate(fl)], fl, np.sqrt(fl), linestyle='', color='black')
        # ax.axhline(mean, color='black', linestyle='--')
        # ax.axhline(regMean, color='red', linestyle='--')

        # ax.set_title('{}\nCCD {}, Ap {} - Poissonian errors - mean flux = {:.3f} counts/s'.format(
        #     obs_T.iso, CCD, '1', mean))
        # ax.set_xlabel('Frame')
        # ax.set_ylabel('Flux, counts/s')
        # plt.show()

        # If we have more than one star, handle that
        if len(ap) > 1:
            for comp in ap[1:]:
                star = data.tseries(CCD, comp)

                # star counts/s
                fl = star.y / data[CCD]['Exptim']

                # Calculate the mean apparent magnitude of the star above the atmosphere
                mag = robust_mag(fl)
                # printer("  Pre-ext correct: CCD {}, Ap {}, mag: {:.3f}".format(CCD, comp, mag))
                mags.append(mag)


        mags = np.array(mags)

        # Add the light lost to atmosphere back in
        printer("  CCD {} extinction: {:.3f} mags".format(CCD, ex*airmass))
        mags = mags - (ex*airmass)


        all_mags[CCD] = np.array(mags)
        del mags
    return all_mags

def get_comparison_magnitudes(std_fname, comp_fname, std_coords, comp_coords,
                                std_mags, obsname, ext):
    '''
    Takes two .log files, one containing the reduction of a standard star and the other the reduction of
    the target frame, using the same settings (aperture size, extraction method, etc.). Uses this to
    compute the apparent magnitudes of comparison stars in comp_fname

    Requires the RA and Dec to correct for airmass.

    Arguments:
    ----------
    std_fname: str
        file containing the standard observations
    comp_fname: str
        file containing the target frame, reduced with identical parameters as the standard
    std_coords: str
        String containing the RA and Dec of the standard, in a style astropy can take.
    comp_coords: str
        String containing the RA and Dec of the target frame, in a style astropy can take.
    std_mags: list
        list containing the SDSS magnitude of the standard in each CCD, inorder
    obsname: str
        Observatory name
    ext: list
        Extinction coefficients, in order of CCD

    Returns:
    --------
    '''
    printer("\n\n--- Extracting comparison star SDSS magnitudes from the file '{}'---".format(comp_fname))
    printer("     using the standard star found in {}\n".format(std_fname))

    std_mags = np.array(std_mags)

    standard_data = hcam.hlog.Hlog.from_ascii(std_fname)
    comp_data     = hcam.hlog.Hlog.from_ascii(comp_fname)

    # Extract the instrumental magnitudes of the standard
    instrumental_std_mags = get_instrumental_mags(standard_data, std_coords, obsname, ext)

    # Convert the dict recieved into an array, so that we have the zero points [r, g, b, ..] in CCD order
    instrumental_std_mags = [instrumental_std_mags[str(i+1)][0] for i, _ in enumerate(instrumental_std_mags)]

    # The zero points are the difference between observed and expected.
    zero_points = instrumental_std_mags - std_mags


    # Get the comparison instrumental mags, in the taget frame
    instrumental_comp_mags = get_instrumental_mags(comp_data, comp_coords, obsname, ext)

    # Discard the variable star magnitude, and do the zero point correction
    apparent_comp_mags = {}
    for i, CCD in enumerate(instrumental_comp_mags):
        ccd = str(i+1)
        instrumental_comp_mags[ccd] = instrumental_comp_mags[ccd][1:]
        apparent_comp_mags[ccd] = instrumental_comp_mags[ccd] - zero_points[i]



    printer("-----------------  STANDARD  -----------------")
    printer("\n  Standard star instrumental magnitudes: ")
    for i, m in enumerate(instrumental_std_mags):
        printer("    CCD {}: {:3.3f}".format(i+1, m))

    printer("\n  Standard Star SDSS magnitudes:")
    for i, m in enumerate(std_mags):
        printer("    CCD {}: {:3.3f}".format(i+1, m))

    printer("\n  Zero points in each band (in order of CCD, will be subtracted from the inst. mags):")
    for i, m in enumerate(zero_points):
        printer("    CCD {}: {:3.3f}".format(i+1, m))

    printer("\n----------------- COMPARISON -----------------")
    printer("\n  Comparison star instrumental magnitudes:")
    for i, _ in enumerate(instrumental_comp_mags):
        CCD = str(i+1)
        printer("    CCD {}: {}".format(CCD,
            np.array2string(instrumental_comp_mags[CCD], precision=3) ))

    printer("\n  Comparison star apparent magnitudes:")
    for i, _ in enumerate(instrumental_comp_mags):
        CCD = str(i+1)
        printer("    CCD {}: {}".format(CCD,
            np.array2string(apparent_comp_mags[CCD], precision=3) ))

    printer('\n  --- Done getting magnitudes ---\n\n')
    return apparent_comp_mags