#!/usr/bin/env python3

# from distutils.core import setup
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
with open("requirements.txt", 'r') as f:
    install_req = f.readlines()

setuptools.setup(
    name='calphot',
    version='1.6',
    description='Set of codes for flux-calibrating *CAM data, when doing relative or absolute photometry',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='James Wild',
    author_email='jwild2@sheffield.ac.uk',
    url='https://github.com/wildjames/cal_phot',
    packages=setuptools.find_packages(),
    package_data={"calphot": ["tab08.dat.txt"]},
    install_requires=install_req,
    python_requires='>=3.6',
    scripts=[
          'scripts/bin_data',
          'scripts/cal_phot',
          'scripts/calc_extinction',
          'scripts/comparison_mags',
          'scripts/splitLogFile'
    ]
)
