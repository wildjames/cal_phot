from ruamel import yaml
import argparse

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


    defaults = {
        'oname': "Reduced_system",
        'get_eclipse_times': False,
        'fit_ephemeris': False,
    }

    with open(yaml_fname) as yaml_file:
        input_dict = yaml.safe_load(yaml_file)

    for key, value in defaults.items():
        if key not in input_dict.keys():
            print("I have no value for {}. Using the default of [{}]".format(key, value))
            input_dict[key] = value

    pprint(input_dict)

    # Scrub the input dict to be in the correct formats.
    # All keys lowercase
    keys = list(input_dict.keys())
    for key in keys:
        input_dict[key.lower()] = input_dict[key]

    # These items must be Boolean
    bool_keys = [
        'sdss',
        'get_eclipse_times',
        'fit_ephemeris',
        'extract_data',
        'flux_calibrate'
    ]
    # Same, for floats
    float_keys = [
        'period',
        't0',
    ]


    #### Do the conversions ####
    ## Bools
    # Top level input dict
    for key in bool_keys:
        if key in input_dict.keys():
            input_dict[key] = bool(input_dict[key])
    # Containers
    outer_keys = list(input_dict.keys())
    for key in input_dict.keys():
        value = input_dict[key]

        if type(value) == dict:
            print("{} is a dict!".format(key))
            inner_keys = list(value.keys())
            for key in bool_keys:
                if key in inner_keys:
                    value[key] = bool(value[key])

    ## Floats
    # Top level input dict
    for key in bool_keys:
        if key in input_dict.keys():
            input_dict[key] = float(input_dict[key])
    # Containers
    outer_keys = list(input_dict.keys())
    for key in input_dict.keys():
        value = input_dict[key]

        if type(value) == dict:
            print("{} is a dict!".format(key))
            inner_keys = list(value.keys())
            for key in bool_keys:
                if key in inner_keys:
                    value[key] = float(value[key])

    pprint(input_dict)


    # Do what the user wants us to
    is_SDSS = input_dict['sdss']

    if is_SDSS:
        pass
    else:
        pass
