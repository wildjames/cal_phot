from ruamel import yaml
import argparse
from pathlib import Path

from getEclipseTimes import getEclipseTimes
from extractData import extract_data
from fitEphem import fitEphem
from plotAll import plot_all

from logger import printer, header


# REMOVEME
from pprint import pprint

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
    }
    for key, value in payload_defaults.items():
        for payload_key, payload in to_extract.items():
            if key not in payload.keys():
                print("{} has no value for {}. Using the default of [{}]".format(payload_key, key, value))
                payload[key] = value

    pprint(input_dict)
    print("\n\n")
    pprint(to_extract)


    # Information gathering
    is_SDSS = input_dict['sdss']
    do_fit_ephem = input_dict['fit_ephemeris']
    target_coords = input_dict['target_coords']
    directory = input_dict['directory']
    T0 = input_dict['t0']
    period = input_dict['period']

    for key, payload in to_extract.items():
        print("Data: {}".format(key))
        pprint(payload)
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

    if do_fit_ephem:
        T0, period = fitEphem(directory, T0, period)
