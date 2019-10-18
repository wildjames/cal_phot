from ruamel import yaml
import argparse
from pathlib import Path
import os

from getEclipseTimes import getEclipseTimes
from extractData import extract_data
from fitEphem import fitEphem
from plotAll import plot_all

from logger import printer, header

#TODO: Write a help file

if __name__ in "__main__":
    parser = argparse.ArgumentParser("YAML input method.")
    parser.add_argument(
        "input_file",
        help="The input YAML file to be computed"
    )
    args = parser.parse_args()
    yaml_fname = args.input_file

    with open(yaml_fname) as yaml_file:
        input_dict = yaml.safe_load(yaml_file)

    # Copy the input file to the Calibration.txt
    header(yaml_file)

    # All keys lowercase
    keys = list(input_dict.keys())
    for key in keys:
        input_dict[key.lower()] = input_dict[key]


    # What data do I have to extract:
    to_extract = input_dict['extract']


    # Set default values
    global_defaults = {
        'directory': '.',
        'fit_ephemeris': False,
    }
    for key, value in global_defaults.items():
        if key not in input_dict.keys():
            input_dict[key] = value

    payload_defaults = {
        'oname': "Reduced_system",
        'get_eclipse_times': False,
        'fit_ephemeris': False,
        'flux_calibrate': True
    }
    for key, value in payload_defaults.items():
        for payload_key, payload in to_extract.items():
            if key not in payload.keys():
                print("{} has no value for {}. Using the default of [{}]".format(payload_key, key, value))
                payload[key] = value

    # Information gathering
    is_SDSS = input_dict['sdss']
    do_fit_ephem = input_dict['fit_ephemeris']
    target_coords = input_dict['target_coords']
    directory = input_dict['directory']
    T0 = input_dict['t0']
    period = input_dict['period']

    # Create the working directory if needed
    if not os.path.isdir(directory):
        os.mkdir(directory)


    # Do the eclipse times, where needed
    for key, payload in to_extract.items():
        print("Data: {}".format(key))
        print("\n\n")

        do_get_ecl_times = payload['get_eclipse_times']

        if do_get_ecl_times:
            print("I want to get eclipse times for {}".format(key))

            obsname = payload['observatory']
            try:
                fnames = payload['logfiles']
            except KeyError:
                print("Searching for log files...")
                globbed = Path('.').glob("**/*.log")
                fnames = list(globbed)
                fnames = [str(f) for f in fnames]
                for fname in fnames:
                    print(fname)

            printer("Getting eclipse times from data...")

            getEclipseTimes(fnames, target_coords, obsname, myLoc=directory)

    # Fit ephemeris
    if do_fit_ephem:
        T0, period = fitEphem(directory, T0, period)

    extracted_files = []
    # Extract the data for each payload
    for key, payload in to_extract.items():
        oname = payload['oname']
        observatory = payload['observatory']
        instrument = payload['inst']

        fnames = payload['logfiles']
        no_calibration = not payload['flux_calibrate']

        if not payload['extract_data']:
            continue

        if no_calibration:
            written_files = extract_data(
                oname, target_coords, observatory, T0, period,
                inst=instrument, SDSS=True, fnames=fnames,
                myLoc=directory
            )
            extracted_files += written_files

        elif is_SDSS:
            written_files = extract_data(
                oname, target_coords, observatory, T0, period,
                inst=instrument, SDSS=is_SDSS, myLoc=directory,
                fnames=fnames
            )
            extracted_files += written_files

        else:
            comparisons = payload['comparison logfiles']
            std_logfile = payload['standard logfile']
            std_coords = payload['standard coords']
            std_mags = payload['standard mags']

            written_files = extract_data(
                oname, target_coords, observatory, T0, period,
                SDSS=is_SDSS, myLoc=directory, fnames=fnames,
                comp_fnames=comparisons, inst=instrument,
                std_fname=std_logfile, std_coords=std_coords, std_mags=std_mags
            )

    if len(written_files):
        print("I created the following files:")
        for f in written_files:
            print(" -> {}".format(f))