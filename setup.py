#!/usr/bin/env python3

from distutils.core import setup

setup(name='cam_calphot',
      version='1.0',
      description='Set of codes for flux-calibrating *CAM data, when doing relative or absolute photometry',
      author='James Wild',
      author_email='jwild2@sheffield.ac.uk',
      url='https://github.com/wildjames/cal_phot',
      packages=['cam_calphot'],
      scripts=['scripts/bin_data', 'scripts/cal_phot', 'scripts/calc_extinction', 'scripts/comparison_mags', 'scripts/splitLogFile'],
     )
