import os
from pathlib import Path

import hipercam as hcam
import matplotlib.pyplot as plt
import numpy as np
from astropy import units
from astropy.coordinates import SkyCoord

from getEclipseTimes import tcorrect

fnames = [str(p) for p in Path(".").glob("**/*.log")]

print("This directory contains these calib files:")
for fname in fnames:
    print(fname)


ephem_fname = "EPHEMERIS/ephemeris_chain.txt"
print("\nI need to know the ephemeris of these data. ")
print("I'll check for {}\n".format(ephem_fname))

if os.path.isfile(ephem_fname):
    # Grab from an ephemeris chain. format is
    # [walker No.] [T0] [P] <...errors...>
    data = np.genfromtxt(ephem_fname, delimiter=' ').T
    T0 = data[1]
    P = data[2]

    # means please
    T0 = np.percentile(T0, 50)
    P = np.percentile(P, 50)
else:
    print("Couldn't find ephemeris data automatically!")
    T0 = input("T0: ")
    P  = input("P:  ")

    T0 = float(T0)
    P = float(P)

print("  T0 = {:.7f}".format(T0))
print("  P  = {:.7f}\n".format(P))

RA = input("Please enter star RA: ")
Dec = input("Please enter star DEC: ")

star = SkyCoord(
    ' '.join([RA,Dec]),
    frame='icrs',
    unit=(units.hour, units.deg)
)



if not os.path.isdir("MCMC_LIGHTCURVES/UNCALIBRATED"):
    os.mkdir("MCMC_LIGHTCURVES/UNCALIBRATED")

fig, ax = plt.subplots()

for fname in fnames:
    print('\nExtracting from {}'.format(fname))

    try:
        data = hcam.hlog.Hlog.read(fname)

        # Check we actually have data.
        if len(data.apnames['1']) < 2:
            print("Failed to extract two apertures from the logfile")
            print("File has apertures: {}".format(data.apnames['1']))
            raise Exception

        target = data.tseries('1', '1') / data.tseries('1', '2')

    # If we fail reading in the log file with the Hipercam software, use the
    # Ultracam version.
    except Exception:
        print("Reading with the ULTRACAM pipeline: {}".format(fname))
        data = hcam.hlog.Hlog.rulog(fname)

        print(data.apnames)
        # Check we actually have data.
        if len(data.apnames['1']) < 2:
            print("Failed to extract two apertures from the logfile")
            print("File has apertures: {}".format(data.apnames['1']))

            print("Skipping file: {}".format(fname))
            continue


    try:
        target = data.tseries('1', '1') / data.tseries('1', '2')
    except Exception as e:
        print(e)
        continue

    observatory = input("Where were these data taken? (site name or lat, lon): ")
    target = tcorrect(target, star, observatory)

    # Phase fold, dumbly
    phi = (target.t - T0) / P
    phi = (phi - 0.5) % 1

    # Scale the flux to means
    flx = target.y / target.ymean()
    err = target.ye / target.ymean()

    ax.errorbar(
        phi, flx, err,
        fmt='-', alpha=0.5,
        label=os.path.split(fname)[-1]
    )

    name = os.path.split(fname)[1].replace('.log', '')
    oname = "MCMC_LIGHTCURVES/UNCALIBRATED/{}.uncalib".format(name)
    with open(oname, 'w') as f:
        f.write("# This is an uncalibrated flux/phase file generated from\n")
        f.write("# {}\n".format(fname))
        for p, fl, fe in zip(phi, flx, err):
            f.write("{} {} {}\n".format(p, fl, fe))

ax.legend()
plt.show()
