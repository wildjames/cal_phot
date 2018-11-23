#!/usr/bin/env python3

import numpy as np
import sys
import os

from getKappa import getKappa
from getEclipseTimes import getEclipseTimes
from combineData import combineData


class Interpreter:
    def __init__(self, inFile=None, prompt=[]):
        # Initialise variables
        self.directory = os.path.curdir
        self.T0 = None
        self.ext = 0.161
        self.binsize = 1
        self.anal_new = True
        self.fnames = None
        self.SDSS = 0

        print("I will use the following working directory: {}\n".format(os.path.curdir+'/Reduced_Data/'))
        if not os.path.exists(os.path.curdir+'/Reduced_Data/'):
            os.mkdir(os.path.curdir+'/Reduced_Data/')


        # parse individual commands
        for line in prompt:
            self.parse(line)
        

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

    def test(self, args):
        print('Entered the testing function, with the arguments:\n  {}'.format(args))
        return None

    def getKappa(self):
        try:
            obsname = self.obsname
            stdLogfile = '/'.join([self.directory, self.stdLogfile])
            mags = self.mags
            coords = self.coords
            ext = self.ext

            print("Computing instrumental magnitude to SDSS magnitude corrections...")
            
            self.kappas = getKappa(stdLogfile, coords, obsname, mags, ext)
            
            print("Computed the following correction factors:\n  r': {:.3f}\n  g': {:.3f}\n  u': {:.3f}\n".format(
                self.kappas[1], self.kappas[2], self.kappas[3]
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
            coords = self.coords
            period = self.period
            obsname = self.obsname
            anal = self.anal_new


            print("Refining eclipse times:\n  Period: {:.3f}\n  Directory: {}\n".format(self.period, self.directory))
            self.T0, self.period = getEclipseTimes(coords, obsname, period, analyse_new=anal, T0=self.T0, myLoc=self.directory)
            print("Using the following ephemeris:\n  T0: {:-6.8f}\n  Period: {:-6.8f}\n".format(
                self.T0, self.period
            ))
        except AttributeError:
            print("I don't have enough data for ephemeris calculation!")
            print('I need the following:')
            print('  - Target coordinates')
            print('  - Initial period guess')
            print('  - Observatory name')
            exit()

    def combineData(self):
        # try:
        oname     = self.oname
        coords    = self.coords
        obsname   = self.obsname
        T0        = self.T0
        period    = self.period
        binsize   = self.binsize
        myLoc     = self.directory
        fnames    = self.fnames
        
        print("Combining, calibrating, and plotting data...")
        if self.SDSS:
            combineData(oname, coords, obsname, T0, period, SDSS=True, binsize=binsize, myLoc=myLoc, fnames=fnames)
        else:
            ref_kappa = self.kappas
            combineData(oname, coords, obsname, T0, period, ref_kappa=ref_kappa, binsize=binsize, myLoc=myLoc, fnames=fnames)
        # except AttributeError:
        #     print("I don't have enough data to do the data processing!")
        #     print("I failed to collect one or more of the following:")
        #     print("  - oname")
        #     print("  - coords")
        #     print("  - ref_kappa")
        #     print("  - observatory name")
        #     print("  - T0")
        #     print("  - period")
        #     print("  - binsize")
        #     print("  - directory")
        #     exit()

    def parse(self, line):
        line = line.strip()

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
            self.obsname = args[0]
            print(f"Observing location: {self.obsname}")
        
        elif command == 'coords':
            # Changes the coordinates of the object you're about to talk about.
            ## TODO: Evaluate -- maybe have separate std and targ coords to ensure users don't forget to change it?
            if len(args) < 2:
                print("I didn't seem to get the right RA and Dec format! Please use:\n  RA - HH:MM:SS.SS\n  DEC - DD:MM:SS.SS\n")
                pass
            else:
                self.coords = '{} {}'.format(args[0], args[1])
                self.coords.replace(':', ' ')
            print("Using the following star coordinates:\n  RA:  {}\n  Dec: {}".format(args[0], args[1]))
        
        elif command == 'directory':
            self.directory = ''.join(args)
            print(f"Working from directory: {self.directory}")

        elif command == 'extinction':
            self.ext = float(args[0])
            print(f"Extinction coefficient: {self.ext}")
        
        elif command == 'plotall':
            self.plotall = True


        # SDSS field observations calibration
        
        elif command == 'sdss':
            # Toggle SDSS field
            self.SDSS = args[0][0]=='1' or args[0][0].lower()=='y'
            print(f"Are we in the SDSS field? [{self.SDSS}]")
        
        elif command == 'sdss_file':
            self.ref_file = args[0]
            print(f"The SDSS reference star RA and Dec are contained in the file '{self.ref_file}'")


        # getKappa stuff
        elif command == 'getkappa':
            self.getKappa()
        
        elif command == 'stdlogfile':
            self.stdLogfile = args[0]
            print("Using the standard star in this log file,\n  {}".format(self.stdLogfile))
        
        elif command == 'stdmags':
            # Must be in the format <r' g' u'>
            if len(args) < 3:
                print("I didn't get enough magnitudes for the standard star!")
            else:
                self.mags = [float(m) for m in args[:3]]
                print("The standard star has the following apparent magnitudes:")
                print("  r': {:.3f}\n  g': {:.3f}\n  u': {:.3f}".format(
                    self.mags[0], self.mags[1], self.mags[2]
                ))
        
        elif command == 'kappa_corr':
            self.kappas = [np.nan, np.nan, np.nan, np.nan]
            self.kappas[1:] = [float(x) for x in args[:3]]
            print("Using the Kappa Corrections:")
            print("  r': {}\n  g': {}\n  u': {}".format(
                self.kappas[1], self.kappas[2], self.kappas[3]
            ))



        # getEclipseTimes
        elif command == 'geteclipsetimes':
            self.getEclipseTimes()
        
        elif command == 'period':
            self.period = float(args[0])
            print(f"Using the period: {self.period}")

        elif command == 't0':
            self.T0 = float(args[0])
            print(f"Using the T0: {self.T0}")

        elif command == 'analyse_new_eclipses':
            self.anal_new = args[0].lower() in ['y', '1']
            if self.anal_new:
                print("I will analyse any data I find for eclipses.")
            else:
                print("Using historical data only for ephemeris fitting.")


        # combineData
        elif command == 'combinedata':
            self.combineData()
        elif command == 'oname':
            self.oname = args[0]
        elif command == 'binsize':
            self.binsize = int(args[0])
        elif command == 'logfiles':
            # Read in logfilenames, terminated by an empty line, i.e. in the format:
            # logfiles
            # file1
            # file2
            # file3 
            #
            # <continue>
            self.fnames = []

            line = self.inFile.readline().strip()
            while line!='':
                self.fnames.append(line)
                line = self.inFile.readline().strip()

            print("Using the following logfiles:")
            for fname in self.fnames:
                print("- {}".format(fname))



        # Unknown command handler
        elif command not in ['', None]: # I don't want to be told about every blank line...
            print("Unknown command!")
            print(f"- {command}")
            print(f"- {args}")


inf = None
prom = []
f = sys.argv

if len(f) == 1:
    interp = Interpreter(prompt=['help'])

if os.path.isfile(f[1]):
    infile = f[1]
    interp = Interpreter(inFile=infile, prompt=prom)
else:
    interp = Interpreter(prompt=f[1:])
