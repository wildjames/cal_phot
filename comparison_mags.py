from constructReference import get_comparison_magnitudes
import hipercam as hcam

if __name__ in "__main__":
    import argparse

    parser = argparse.ArgumentParser("Calculate the apparent magnitude of each aperture in a run, calibrated from a standard star observation")
    parser.add_argument(
        "std_fname",
        help="The logfile reduction of the standard star"
    )
    parser.add_argument(
        "comp_fname",
        help="The logfile reduction of the comparison stars, with the same settings as the standard"
    )
    parser.add_argument(
        "obsname",
        help="The name of the observing site"
    )

    args = parser.parse_args()
    std_fname = args.std_fname
    comp_fname = args.comp_fname
    obsname = args.obsname

    std_coords = input("Please enter the standard RA, Dec: ")
    comp_coords = input("Please enter the comparison stars' RA, Dec (for airmass calc): ")

    std_data = hcam.hlog.Hlog.read(std_fname)
    aps = std_data.apnames
    nCCD = len(aps.keys())
    print("The standard reduction in {} has {} CCDs".format(std_fname, nCCD))

    std_mags = []
    print("\nPlease enter standard mags, in CCD order:")
    while len(std_mags) != nCCD:
        print("Currently, I have mags: {}".format(std_mags))
        mag = input("Enter a mag: ")
        mag = float(mag)
        std_mags.append(mag)
        print()

    ext = []
    print("\n\nPlease give me the extinction of each CCD, in order:")
    while len(ext) != nCCD:
        print("Currently, I have extinctions (mags/airmass): {}".format(ext))
        ex = input("Enter k_(ext): ")
        ex = float(ex)
        ext.append(ex)

    mags = get_comparison_magnitudes(
        std_fname, comp_fname,
        std_coords, comp_coords,
        std_mags,
        obsname,
        ext
    )

    print("Here's the bloc that you'll want to put in the config YAML for cal_phot:")
    for i, CCD in mags.items():
        string = ["{:.3f}".format(mag) for mag in CCD]
        print("  - [{}]".format(','.join(string)))