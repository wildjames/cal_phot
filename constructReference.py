import json
import requests
from pprint import pprint
import os

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

    print("Getting reference SDSS magnitudes from '{:s}'".format(fetchFname.split('/')[-1]))

    CCDs = ['1', '2', '3']

    fetchme = {
        "1": [],
        "2": [],
        "3": []
    }
    
    if not os.path.isfile(fetchFname):
        # Check that we;ve not been given a directory with a 'coord_list' file in it:
        test = fetchFname.split('/')
        test.append('coord_list.coords')
        test = '/'.join(test)
        if os.path.isfile(test):
            print("I found the following file in that directory:")
            print(test)
            fetchFname = test
        # If not, then we have no file. Create a template.
        else:
            with open(fetchFname, 'w') as f:
                f.write('''<CCD 1 reference 1 RA> <CCD 1 reference 1 Dec>
<CCD 1 reference 2 RA> <CCD 1 reference 2 Dec>

<CCD 2 reference 1 RA> <CCD 2 reference 1 Dec>
<CCD 2 reference 2 RA> <CCD 2 reference 2 Dec>
<CCD 2 reference 3 RA> <CCD 2 reference 3 Dec>

<CCD 3 reference 1 RA> <CCD 3 reference 1 Dec>
<CCD 3 reference 2 RA> <CCD 3 reference 2 Dec>''')
            print("Couldn't find that co-ordinate list! I created a template for you at the location you gave me.")
            return None

    with open(fetchFname) as f:
        x = 1
        for line in f:
            if x > 3:
                break
            if line[0] == '#':
                # print(line.strip())
                pass
            elif len(line.split()) != 2:
                x += 1
            else:
                fetchme[str(x)].append(line.split())

    toWrite = {
        '1':{},
        '2':{},
        '3':{}
    }

    for CCD in CCDs:
        # print('CCD {}'.format(CCD))
        # Grab the list of coordinates we want to query
        coords = fetchme[CCD]

        for i, coord in enumerate(coords):
            ra, dec = coord

            # print('  ra: {}, dec: {}'.format(ra, dec))

            # Construct the URL we're gonna post. First define what DB we want to search
            url  = 'http://skyserver.sdss.org/dr14/SkyserverWS/SearchTools/RadialSearch?'
            url += 'ra={}&dec={}&'.format(ra, dec)
            # I'm using a radial search, this is the readius of that search
            url += 'radius=0.1&'
            # Which coord system are we using. Can't imagine I'll want galactic...
            url += 'whichway=equitorial&'
            # If I'm getting more than 5 results, I've probably picked a crappy reference anyway...
            url += 'limit=5&'
            url += 'format=json&'
            url += 'whichquery=imaging'

            resp = requests.post(url)

            # pprint(resp.json())

            results = resp.json()[0]['Rows']
            if len(results) >= 5:
                print('''You got a lot of results from the query. 
You probably want to go and manually add this star to the JSON!''')
            if len(results) > 1:
                print('''
--------------------------------------------
More than one object found at that location!'''
                )
                # Get the user to pick one:
                for m, line in enumerate(results):
                    print("{}\n  RA: {}, Dec: {}\n  u: {}\n  g: {}\n  r: {}".format(
                        m, line['ra'], line['dec'], 
                        line['u'], line['g'], line['r']
                        )
                    )
                n = input("Which object to use?: ")
                target = results[int(n)]
                print('--------------------------------------------')
            elif len(results) == 1:
                target = results[0]
                ra = deg2arcsec(target['ra'], ra=True)
                dec = deg2arcsec(target['dec'], ra=False)
                # print("    Found one result:\n      ra: {}, dec: {}\n        u: {}\n        g: {}\n        r: {}".format(
                #         ra, dec, 
                #         target['u'], target['g'], target['r']
                #         )
                #     )
            else:
                print('''ERROR! Found no targets at this location!!)
Try broadening the search radius in this script,
and make sure that your targets are definitely in the SDSS field!''')
                continue

            # pprint(target)
            toWrite[CCD][str(i+2)] = target

    print("Done!\n")
    return toWrite