import json
import os

import bs4 as bs
import hipercam as hcam
import numpy as np
import requests
from astropy import coordinates as coord
from astropy import time
from astropy import units as u
from astropy.coordinates import AltAz
from astropy.stats import sigma_clipped_stats

try:
    from .logger import printer
except ImportError:
    def printer(string, end='\n'):
        print(string, end=end)




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
    '''Converts a list of fluxes to magnitudes'''
    try:
        mean, median, sigma = sigma_clipped_stats(cps, maxiters=2, sigma=3)
    except:
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

def convert_kg5(sdss_result):
    r = sdss_result['r']
    g = sdss_result['g']

    KG5 = g  - 0.2240*(g-r)**2 - 0.3590*(g-r) + 0.0460
    printer("  ** Using the recipe: KG5 = g - 0.2240*(g-r)**2 - 0.3590*(g-r) + 0.0460")
    printer("  ** Computed the KG5 magnitude: {:.3f}".format(KG5))

    return KG5


def construct_reference(fetchFname):
    '''
    Queries the SDSS database for the magnitudes of the stars contained in <fetchFname>, returns them as a dict of lists.

    mags = {
        '1': [ap1, ap2, ap3],      # r' band
        '2': [ap1, ap2],           # g' band
        '3': [ap1, ap2, ap3, ap4]  # u' band
    }

    <fetchFname> is formatted as follows:
        <CCD1 filter> <CCD2 filter> <CCD3 filter>

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

    radius = '0.1' # search radius, arcseconds


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
                f.write("#For the example of a 3-band instrument, use the following format. Similar for N bands")
                f.write("<CCD1 SDSS band> <CCD2 SDSS band> <CCD3 SDSS band>\n")
                f.write("\n")
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
        line = f.readline()
        while line[0] == '#':
            line = f.readline()
        line = line.strip().split()
        bands = ['']
        for b in line:
            bands.append(b)
        f.readline()
        printer('Querying SDSS for the following bands: {}'.format(', '.join(bands[1:])))

        # Construct a dict to read the desired coordinates to
        fetchme = {str(i+1):[] for i, b in enumerate(bands[1:])}
        # Make a list of the CCDs, starting from 1
        CCDs    = [str(i+1)    for i, b in enumerate(bands[1:])]

        x = 1
        for line in f:
            if line[0] == '#':
                # print(line.strip())
                pass
            # If the line doesn't have exactly two things, it's not something I care about.
            elif len(line.split()) != 2:
                x += 1
            else:
                fetchme[str(x)].append(line.split())

    printer("  Looking for the following the coordinates:")
    printer("  {:^15s} | {:^15s}".format('RA', 'DEC'))
    for CCD in CCDs:
        coords = fetchme[CCD]
        printer("    {:<13s} | {:<15s}".format('CCD {}'.format(CCD), ''))
        for coord in coords:
            ra, dec = coord
            printer("  {:>15s} | {:<15s}".format(ra, dec))

    toWrite = {}

    for CCD in CCDs:
        printer('-> CCD {}'.format(CCD))
        # Grab the list of coordinates we want to query
        coords = fetchme[CCD]
        band = bands[int(CCD)]

        for index, coord in enumerate(coords):
            ra, dec = coord

            printer('    Searching -> RA, Dec: {}, {} for {} band mag'.format(ra, dec, band))

            # Construct the URL we're gonna post. First define what DB we want to search
            url  = 'http://skyserver.sdss.org/dr14/SkyserverWS/SearchTools/RadialSearch?'
            url += 'ra={}&dec={}&'.format(ra, dec)
            # I'm using a radial search, this is the readius of that search
            url += 'radius={}&'.format(radius)
            # Which coord system are we using. Can't imagine that I'll want galactic
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
                printer("If you want to look in more detail, here's my search URL:\n{}".format(url.replace('json', 'html')))
                # Get the user to pick one:
                for m, line in enumerate(results):
                    printer("{}\n  RA, Dec: {} {}\n  u: {}\n  g: {}\n  r: {}\n  i: {}\n  z: {}".format(
                        m, line['ra'], line['dec'],
                        line['u'], line['g'], line['r'], line['i'], line['z']
                        )
                    )
                n = input("Which object to use?: ")
                while n == '':
                    n = input("Which object to use?: ")
                    if int(n) >= len(results):
                        n = ''
                printer("Chose object {}".format(n), terminal=False)
                target = results[int(n)]
                printer('--------------------------------------------')
            elif len(results) == 1:
                target = results[0]
                ra = deg2arcsec(target['ra'], ra=True)
                dec = deg2arcsec(target['dec'], ra=False)
                printer("    Found one result:\n      ra: {}, dec: {}\n       u: {}\n       g: {}\n       r: {}\n       i: {}\n       z: {}".format(
                        ra, dec,
                        target['u'], target['g'], target['r'], target['i'], target['z']
                        )
                    )

            else:
                printer('ERROR! Found no targets at the location: RA: {}, Dec: {}'.format(target['ra'], target['dec']))
                printer('Try broadening the search radius in this script (was {}),'.format(radius))
                printer('and make sure that your targets are definitely in the SDSS field!')
                raise LookupError

            # append the magnitudes found in [bands] to the output dict.
            if band.lower() == 'kg5':
                target[band] = convert_kg5(target)


            ### I need to check the flags here. Mini webscraper:
            target_entry = 'http://skyserver.sdss.org/DR14/en/tools/explore/summary.aspx?id={}'.format(target['objid'])
            printer("    Here's the entry on the skyserver:\n    {}".format(target_entry))
            # Lovely soup
            resp = requests.post(target_entry)
            soup = bs.BeautifulSoup(resp.text, features='lxml')

            # The flags are stored in a table, with the first column having the string "Flags"...
            tables = soup.find_all('table')

            foundFlags = False
            for table in tables:
                # Scan through for the table that keeps the flags. USE HTML HEADERS PEOPLE
                if 'Flags' in table.text:
                    # Yay!
                    foundFlags = True
                    FLAGS = table.text.strip().split()[1:]
                    if len(FLAGS):
                        printer("    This star has the following flags:")
                        for flag in FLAGS:
                            printer("      -> {}".format(flag))

                        if 'SATURATED' in FLAGS:
                            printer("THIS STAR HAS SATURATED SDSS, AND WILL NOT GIVE AN ACCURATE FLUX.")
                            printer("I'll infer its magnitude from other comparisons..")
                            target[band] = np.nan

                        stop = input("Hit enter to continue if these are okay, 'q' to stop the script: ") + ' '
                        if 'q' in stop.lower():
                            printer("User terminated during SDSS star lookups")
                            quit()
                        break

            if not foundFlags:
                # Boooooo
                printer("Didn't find the flags table... Check this one manually!")
                printer("\n")
                input("Cont...")

            printer("    -> Using the {} band magnitude for star {}: {:.3f}".format(band, index, target[band]))
            printer("\n\n")

            try:
                toWrite[CCD].append(
                    target[band]
                )
            except KeyError:
                toWrite[CCD] = [ target[band] ]

        toWrite[CCD] = np.array(toWrite[CCD])

    printer("Got all reference stars for this file!\n")
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

    ext: iterable
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
        try:
            observatory = coord.EarthLocation.of_site(obsname)
            printer("Observatory successfully retrieved from site name")
        except:
            obsname = obsname.split(',')
            if len(obsname) != 2:
                printer("  The (lat, lon) MUST!!! be comma separated!")
                exit()
            lat, lon = obsname
            printer("  Earth location from latitude, longitude: {}, {}".format(lat, lon))
            observatory = coord.EarthLocation.from_geodetic(lat=lat, lon=lon)

        star_loc = coord.SkyCoord(
            coords,
            unit=(u.hourangle, u.deg)
        )

        # I want altitude converted to zenith angle. Airmass is roughly constant over
        # a single eclipse so only do it once to save time.
        obs_T = data.tseries('1', '1').t
        obs_T = time.Time(obs_T, format='mjd')

        # Define the altAz frame, containing the time and location
        altAz_frame = AltAz(obstime=obs_T, location=observatory)
        star_loc_AltAz = star_loc.transform_to(altAz_frame)

         # Compute the airmass, at the time of the first frame
        zenith_angle = 90 - star_loc_AltAz.alt.deg
        zenith_angle_rad = np.deg2rad(zenith_angle)
        airmass = 1. / np.cos(zenith_angle_rad)

        printer(
            "  For the observations starting at {} and ending at {}...".format(
                obs_T[0].iso, obs_T[-1].iso
            )
        )
        printer(
            "   -> Zenith angle starts at {:.3f}, and ends at {:.3f}".format(
                zenith_angle[0], zenith_angle[-1]
            )
        )
        printer(
            "   -> Airmass starts at {:.3f}, ends at {:.3f}".format(
                airmass[0], airmass[-1]
            )
        )
    else:
        printer("  No coordinates or observatory provided, setting airmass to 0.0")
        airmass = [0.0 for _ in data.tseries('1', '1').t]

    ##TODO: The mean airmass is used for now.
    airmass = np.mean(airmass)
    print("Mean airmass: {:.3f}".format(airmass))
    if airmass <= 0:
        printer("Airmass is negative!! We have a problem there!")
        printer("EarthLocation (constructed from {}):".format(obsname))
        printer(str(observatory))
        printer(str(observatory.lat), str(observatory.lon))
        printer("Star location:")
        printer(str(star_loc))
        input("Hit enter to continue... ")

    printer("Getting the INSTRUMENTAL (electron flux) magnitudes for the log file")

    all_mags = {}
    aps = data.apnames
    CCDs = [str(i+1) for i, key in enumerate(aps)]

    if ext is None:
        ext = [0.0 for i in CCDs]
    ext = np.array(ext)

    for CCD in CCDs:
        printer("\n---> Doing CCD {} <---".format(CCD))
        # Get this frame's apertures
        ap = sorted(aps[CCD])

        ex = ext[int(CCD)-1]

        star = data.tseries(CCD, '1')

        exptime = data[CCD]['Exptim']

        # star counts/s
        fl = star.y / exptime

        printer("The first aperture had a mean counts per frame of {:.2f}".format(np.mean(star.y)))
        printer("  and a mean exposure time of {:.3f}".format(np.mean(exptime)))

        # Calculate the mean apparent magnitude of the star above the atmosphere
        mag = robust_mag(fl)
        # star magnitudes
        mags = [mag]

        # If we have more than one star, handle that
        if len(ap) > 1:
            for comp in ap[1:]:
                star = data.tseries(CCD, comp)

                # star counts/s
                fl = star.y / exptime

                # Filter out bad data
                if np.any(star.mask):
                    print("This data has bad flags!")

                printer("Aperture {} had a mean counts per frame of {:.2f}".format(comp, np.mean(star.y)))
                printer("  and a mean exposure time of {:.3f}".format(np.mean(exptime)))

                # Calculate the mean apparent magnitude of the star above the atmosphere
                mag = robust_mag(fl)
                printer("  Pre-ext correct: CCD {}, Ap {}, mag: {:.3f}".format(CCD, comp, mag))
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

    standard_data = hcam.hlog.Hlog.read(std_fname)
    comp_data     = hcam.hlog.Hlog.read(comp_fname)

    # Extract the instrumental magnitudes of the standard
    printer("   -> Computing the instrumental magnitudes of the standard")
    instrumental_std_mags = get_instrumental_mags(standard_data, std_coords, obsname, ext)

    # Convert the dict recieved into an array, so that we have the zero points [r, g, b, ..] in CCD order
    instrumental_std_mags = [instrumental_std_mags[str(i+1)][0] for i, _ in enumerate(instrumental_std_mags)]
    instrumental_std_mags = np.array(instrumental_std_mags)

    # The zero points are the difference between observed and expected.
    zero_points = instrumental_std_mags - std_mags


    # Get the comparison instrumental mags, in the taget frame
    printer("\n\n   -> Computing the instrumental magnitudes of the target frame")
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
