#!/usr/bin/env python3

import numpy as np
import sys
import os

from getKappa import getKappa
from getEclipseTimes import getEclipseTimes
from combineData import combineData


class Interpreter:
    def __init__(self, inFile=None, prompt=[]):
        # Initialise variables. Store args in a dict.
        self.params = {
            'T0': None,
            'directory': os.path.curdir,
            'ext': 0.161,
            'binsize': 1,
            'anal_new': True,
            'fnames': None,
            'SDSS': 0
        }

        # parse prompt commands
        for line in prompt:
            self.parse(line)
        
        # parse file commands
        if inFile != None:
            self.inFile = open(inFile, 'r')
            line = self.inFile.readline()
            self.parse(line)
            while line:
                line = self.inFile.readline()
                self.parse(line)
            
            self.inFile.close()
    
    def help(self):
        print('''
This is an interpreter script that sits on top of a few data HiPERCAM/ULTRACAM reduction scripts I wrote to do calibrated photometry.
As a warning, I keep my data in this kind of structure:
- Object name
-- Observing night 1
-- Observing night 2
If you want to ensure the script works properly, it's best to follow this convention.

Generally we want to follow these steps:
- Are we in the SDSS field? 
-- If we are in the SDSS, I'll perform a database query for your standards, but their coordinates must be supplied to let me do this.
-- If we aren't in the SDSS, you must reduce data for a standard star, and pass the information on where to find that to me.
--- Then, I need the RA and Dec, as well as the SDSS magnitudes of the star and the location it was observed from to get the corrections.
- I need the ephemeris data on the system, T0 and period, in order to fold it.
-- If you have a set of eclipse times, or some eclipse lightcurves, or both, I have a script to refine your period estimate.
- Then, show me where all your logfiles are or let me search for them, and I'll collect and calibrate your lightcurves from counts to 
   fluxes, fold them about the period, and plot you a lightcurve. The reduced data are also written to a .calib file for future use.


---- COMMANDS ----

-- Housekeeping commands --
  test <args>
    Testing command, literally just prints the arguments.
  help
    Print this text
  stop/exit/quit
    Stops the processing
  observatory <name>
    Sets the observing location to the argument it's given. run astropy.EarthLocation.get_site_names() to get a list of valid sites.
  coords <RA> <Dec
    Sets the sky location. This has to be changed if your standard is in a different place than your target!
  directory <dir>
    If we are using sub-folders, this tells the next script where to look for data
  extinction <ext>
    Sets extinction coefficient, in mags/airmass
  writeparams <paramName>
    Writes parameters to file. You probably want this every time! If no paramname is given, writes to 'reductino_params.txt'.

-- Kappa corrections commands --
  stdLogfile <file>
    Tells the kappa correction script where to look for a .log file
  stdMags <r' mag> <g' mag> <u' mag>
    SDSS magnitudes for the standard star
  kappa_corr <r' correction> <g' correction> <u' correction>
    Manually define the correction factors

-- Ephemeris commands --
  period <P>
    Sets the period of the system for folding. If we are doing getEclipseTimes, this forms the initial guess
  T0 <T0>
    Sets the T0 for folding. 

-- Data combination commands --
  oname <String>
    Define the names prefix to use when writing out reduced data
  binsize <int>
    Define the binning after we've folded our data. <1> sets to unbinned.

-- 'Trigger' commands
  getKappa
    Calls the script to calculate the photometric -> magnitude corrections.
    Requires:
     - An hcam logfile containing the photometry of the standard star.
     - Std. star coordinates
     - Std. star SDSS magnitudes
     - An extinction coefficient, in mags/airmass. This defaults to 0.161
     - The name of the observatory where the data were taken. Used to lookup in astropy
  
  getEclipseTimes
    Calls a script that searches the given directory for logfiles, and with the user's help determines 
    the eclipse times of the lightcurves. These are then used to compute the best period for the data.
    Requires: 
     - Target coordiates
     - Initial period guess
     - [OPTIONAL] Previously found value of T0

  combineData
    Searches for logfiles, callibrates the photometry, folds and bins the data.
    Requires: 
     - Target coordinates
     - Kappa corrections
     - Observatory name
     - Ephemeris data, T0 and period
     - Binsize
     - Directory
     - Name template to use while writing out to files

''')

    def get_param(self, pname):
        try:
            p = self.params[pname]
        except KeyError:
            p = None
            print("I couldn't retrieve the parameter {}!".format(pname))
            raise AttributeError
        return p

    def test(self, args):
        print('Entered the testing function, with the arguments:\n  {}'.format(args))
        return None

    def getKappa(self):
        try:
            obsname = self.get_param('obsname')
            mags   = self.get_param('mags')
            coords = self.get_param('coords')
            ext    = self.get_param('ext')
            stdLogfile = '/'.join(
                [self.get_param('directory'), self.get_param('stdLogfile')]
            )
            
            print("Computing instrumental magnitude to SDSS magnitude corrections...")
            
            kappas = getKappa(stdLogfile, coords, obsname, mags, ext)
            self.params['kappas'] = kappas

            print("Computed the following correction factors:\n  r': {:.3f}\n  g': {:.3f}\n  u': {:.3f}\n".format(
                kappas[1], kappas[2], kappas[3]
                )
            )
        except AttributeError:
            print("I don't have enough data to determine Kappa correction!")
            print("I need the following:")
            print("  - Std. star photometry, in hipercam logfile.")
            print("  - Std. star coordinates")
            print("  - Std. star SDSS magnitudes, r', g', u'.")
            print("  - Extinction coefficient")
            print("  - Observatory name")
            exit()

    def getEclipseTimes(self):
        try:
            coords    = self.get_param('coords')
            period    = self.get_param('period')
            obsname   = self.get_param('obsname')
            anal      = self.get_param('anal_new')
            directory = self.get_param('directory')
            T0        = self.get_param('T0')


            print("Refining eclipse times:\n  Period: {:.3f}\n  Directory: {}\n".format(period, directory))

            T0, period = getEclipseTimes(coords, obsname, period, analyse_new=anal, T0=T0, myLoc=directory)
            self.params['T0']     = T0
            self.params['period'] = period
            
            print("Using the following ephemeris:\n  T0: {:-6.8f}\n  Period: {:-6.8f}\n".format(
                T0, period
            ))
        except AttributeError:
            print("I don't have enough data for ephemeris calculation!")
            print('I need the following:')
            print('  - Target coordinates')
            print('  - Initial period guess')
            print('  - Observatory name')
            exit()

    def combineData(self):
        try:
            oname     = self.get_param('oname')
            coords    = self.get_param('coords')
            obsname   = self.get_param('obsname')
            T0        = self.get_param('T0')
            period    = self.get_param('period')
            binsize   = self.get_param('binsize')
            myLoc     = self.get_param('directory')
            fnames    = self.get_param('fnames')
            SDSS      = self.get_param('SDSS')
            
            print("Combining, calibrating, and plotting data...")
            if SDSS:
                combineData(oname, coords, obsname, T0, period, SDSS=True, binsize=binsize, myLoc=myLoc, fnames=fnames)
            else:
                ref_kappa = self.get_param('kappas')
                combineData(oname, coords, obsname, T0, period, ref_kappa=ref_kappa, SDSS=False, binsize=binsize, myLoc=myLoc, fnames=fnames)
        except AttributeError:
            print("I don't have enough data to do the data processing!")
            print("I failed to collect one or more of the following:")
            print("  - oname")
            print("  - coords")
            print("  - ref_kappa")
            print("  - observatory name")
            print("  - T0")
            print("  - period")
            print("  - binsize")
            print("  - directory")
            exit()
        else:
            print("combineData command had a weird error...")

    def parse(self, line):
        line = line.strip()

        # Clean up input for parsing
        if line == '':
            command = None
            args = None
        elif line[0] == '#':
            command = None
            args = None
        elif '#' in line:
            cut = line.index('#')
            line = line[:cut]
            line = line.split(' ')
            command = line[0].lower()
            if len(line) > 1:
                args = [x for x in line[1:] if x not in ['', ',', ' ']]
            else:
                args = None
        else:
            # Split the line into command and arguments, space separated
            line = line.split(' ')
            command = line[0].lower()
            if len(line) > 1:
                args = [x for x in line[1:] if x not in ['', ',', ' ']]
            else:
                args = None

        # print("  Command: {}\n  args: {}".format(command, args))
        
        ## Housekeeping commands
        if   command == 'test':
            self.test(args)
        elif command == 'help':
            self.help()
            exit()
        elif command in ['stop', 'exit', 'quit']:
            print("Stopping.")
            exit()

        ## Actual commands
        # ''Global'' variables
        elif command == 'observatory':
            # Changes the observing location
            obsname = args[0]
            self.params['obsname'] = obsname
            print("Observing location: {}".format(obsname))
        
        elif command == 'coords':
            # Changes the coordinates of the object you're about to talk about.
            ## TODO: Evaluate -- maybe have separate std and targ coords to ensure users don't forget to change it?
            if len(args) < 2:
                print("I didn't seem to get the right RA and Dec format! Please use:\n  RA - HH:MM:SS.SS\n  DEC - DD:MM:SS.SS\n")
                pass
            else:
                coords = '{} {}'.format(args[0], args[1])
                coords.replace(':', ' ')
                self.params['coords'] = coords
                print("Using the following star coordinates:\n  RA:  {}\n  Dec: {}".format(args[0], args[1]))
        
        elif command == 'directory':
            directory = ''.join(args)
            self.params['directory'] = directory
            print("Working from directory: {}".format(directory))

        elif command == 'extinction':
            ext = float(args[0])
            self.params['ext'] = ext
            print("Extinction coefficient: {}".format(ext))
        
        elif command == 'plotall':
            self.params['plotall'] = args[0] in ['y', '1', 'yes', 'true']

        elif command == 'writeparams':
            if args == []:
                paramname = 'reduction_params.txt'
            else:
                paramname = args[0]
            with open(paramname, 'w') as f:
                for key, item in enumerate(self.params):
                    f.write("{} {}\n".format(item, self.params[item]))
            print("Wrote parameters to 'reduction_params.txt'!")


        # SDSS field observations calibration
        
        elif command == 'sdss':
            # Toggle SDSS field
            SDSS = args[0] in ['y', '1', 'yes', 'true']
            self.params['SDSS'] = SDSS
            print("Are we in the SDSS field? [{}]".format(SDSS))
        
        elif command == 'sdss_file':
            ref_file = args[0]
            self.params['ref_file'] = ref_file
            print("The SDSS reference star RA and Dec are contained in the file '{}'".format(ref_file))


        # getKappa stuff
        elif command == 'getkappa':
            self.getKappa()
        
        elif command == 'stdlogfile':
            stdLogfile = args[0]
            self.params['stdLogfile'] = stdLogfile
            print("Using the standard star in this log file,\n  {}".format(stdLogfile))
        
        elif command == 'stdmags':
            # Must be in the format <r' g' u'>
            if len(args) < 3:
                print("I didn't get enough magnitudes for the standard star!")
            else:
                mags = [float(m) for m in args[:3]]
                self.params['mags'] = mags

                print("The standard star has the following apparent magnitudes:")
                print("  r': {:.3f}\n  g': {:.3f}\n  u': {:.3f}".format(
                    mags[0], mags[1], mags[2]
                ))
        
        elif command == 'kappa_corr':
            kappas = [np.nan, np.nan, np.nan, np.nan]
            kappas[1:] = [float(x) for x in args[:3]]
            self.params['kappas']

            print("Using the Kappa Corrections:")
            print("  r': {}\n  g': {}\n  u': {}".format(
                kappas[1], kappas[2], kappas[3]
            ))



        # getEclipseTimes
        elif command == 'geteclipsetimes':
            self.getEclipseTimes()
        
        elif command == 'period':
            period = float(args[0])
            self.params['period'] = period
            print("Using the period: {}".format(period))

        elif command == 't0':
            T0 = float(args[0])
            self.params['T0'] = T0
            print("Using the T0: {}".format(T0))

        elif command == 'analyse_new_eclipses':
            anal_new = args[0].lower() in ['y', '1']
            self.params['anal_new'] = anal_new
            if anal_new:
                print("I will analyse any data I find for eclipses.")
            else:
                print("Using historical data only for ephemeris fitting.")


        # combineData
        elif command == 'combinedata':
            self.combineData()
        elif command == 'oname':
            oname = args[0]
            self.params['oname'] = oname
            print("Using the following filename: {}".format(oname))
        elif command == 'binsize':
            binsize = int(args[0])
            self.params['binsize'] = binsize
            print("Binning data by {}".format(binsize))
        elif command == 'logfiles':
            # Read in logfilenames, terminated by an empty line, i.e. in the format:
            # logfiles
            # file1
            # file2
            # file3 
            #
            # <continue>
            fnames = []

            line = self.inFile.readline().strip()
            while line!='':
                fnames.append(line)
                line = self.inFile.readline().strip()
            self.params['fnames'] = fnames

            print("Using the following logfiles:")
            for fname in fnames:
                print("- {}".format(fname))



        # Unknown command handler
        elif command not in ['', None]: # I don't want to be told about every blank line...
            print("Unknown command!")
            print("- {}".format(command))
            print("- {}".format(args))


inf = None
prom = []
f = sys.argv

if len(f) == 1:
    interp = Interpreter(prompt=['help'])

if os.path.isfile(f[1]):
    infile = f[1]
    interp = Interpreter(inFile=infile)
else:
    interp = Interpreter(prompt=f[1:])
