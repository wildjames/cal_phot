#!/usr/bin/env python3

import numpy as np
import sys
import os

from getKappa import getKappa
from getEclipseTimes import getEclipseTimes
from combineData import combineData
from fitEphem import fitEphem
from plotAll import plot_all

from logger import printer, header

class Interpreter:
    def __init__(self, inFile=None, prompt=False):
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
        self.written_files = []
        
        # parse file commands
        if inFile != None:
            header(inFile)

            self.inFile = open(inFile, 'r')
            line = self.inFile.readline()
            self.parse(line)
            while line:
                try:
                    line = self.inFile.readline()
                    self.parse(line)
                except ValueError:
                    break
            self.inFile.close()
        
        if prompt:
            while True:
                line = input("> ")
                self.parse(line)
    
    def help(self):
        print('''
        This is an interpreter script that sits on top of a few data HiPERCAM/ULTRACAM reduction scripts I wrote to do calibrated photometry.
        As a warning, I keep my data in this kind of structure:
        - Object name
        -- Observing night 1
        -- Observing night 2
        If you want to ensure the script works properly, it's best to follow this convention.
        ''')

    def get_param(self, pname):
        try:
            p = self.params[pname]
        except KeyError:
            p = None
            printer("I couldn't retrieve the parameter {}!".format(pname))
            raise AttributeError
        return p

    def test(self, args):
        printer('Entered the testing function, with the arguments:\n  {}'.format(args))
        return None

    def getKappa(self):
        obsname    = self.get_param('obsname')
        mags       = self.get_param('mags')
        coords     = self.get_param('coords')
        ext        = self.get_param('ext')
        stdLogfile = self.get_param('stdLogfile')
        
        printer("Computing instrumental magnitude to SDSS magnitude corrections...")
        
        kappas = getKappa(stdLogfile, coords, obsname, mags, ext)
        self.params['kappas'] = kappas

        printer("Computed the following correction factors:\n  r': {:.3f}\n  g': {:.3f}\n  u': {:.3f}\n".format(
            kappas[1], kappas[2], kappas[3]
            )
        )
        
    def getEclipseTimes(self):
        coords    = self.get_param('coords')
        obsname   = self.get_param('obsname')
        directory = self.get_param('directory')

        printer("Getting eclipse times from data...")

        getEclipseTimes(coords, obsname, myLoc=directory)
        
    
    def fitEphem(self):
        directory = self.get_param('directory')
        T0        = self.get_param('T0')
        period    = self.get_param('period')

        printer("Fitting eclipse times to refine period and T0 parameters")

        T0, period = fitEphem(directory, T0, period)

    def combineData(self):
        oname     = self.get_param('oname')
        coords    = self.get_param('coords')
        obsname   = self.get_param('obsname')
        T0        = self.get_param('T0')
        period    = self.get_param('period')
        binsize   = self.get_param('binsize')
        myLoc     = self.get_param('directory')
        fnames    = self.get_param('fnames')
        SDSS      = self.get_param('SDSS')
        
        printer("Combining, calibrating, and plotting data...")
        if SDSS:
            written_files = combineData(oname, coords, obsname, T0, period, SDSS=True, binsize=binsize,
            myLoc=myLoc, fnames=fnames)
        else:
            # Retrieve the SDSS-matching reductions for each night
            comparisons = self.get_param('comparisonfnames')
            stdLogfile  = self.get_param('stdLogfile')
            stdCoords   = self.get_param('stdcoords')
            stdMags     = self.get_param('mags')
            
            # ref_kappa = self.get_param('kappas')
            written_files = combineData(oname, coords, obsname, T0, period, SDSS=False,
                binsize=binsize, myLoc=myLoc, fnames=fnames, comp_fnames=comparisons, 
                std_fname=stdLogfile, std_coords=stdCoords, std_mags=stdMags)

        self.written_files += written_files
        
        printer("So far, I've written the following files:")
        for f in self.written_files:
            printer("-> {}".format(f))

    def parse(self, line):
        line = line.strip()

        # Clean up input for parsing
        if line == '':
            command = None
            args = []
        elif line[0] == '#':
            command = None
            args = []
        elif '#' in line:
            cut = line.index('#')
            line = line[:cut]
            line = line.split(' ')
            command = line[0].lower()
            if len(line) > 1:
                args = [x for x in line[1:] if x not in ['', ',', ' ']]
            else:
                args = []
        else:
            # Split the line into command and arguments, space separated
            line = line.split(' ')
            command = line[0].lower()
            if len(line) > 1:
                args = [x for x in line[1:] if x not in ['', ',', ' ']]
            else:
                args = []

        # print("  Command: {}\n  args: {}".format(command, args))
        
        ## Housekeeping commands
        if   command == 'test':
            self.test(args)
        elif command == 'help':
            self.help()
            exit()
        elif command in ['stop', 'exit', 'quit']:
            printer("Stopping.")
            print("\n\n\n\n")
            try:
                self.inFile.close()
            except AttributeError:
                exit()

        ## Actual commands
        # '''Global''' variables
        elif command == 'directory':
            directory = args[0]
            self.params['directory'] = directory
            printer("Working from directory: {}".format(directory))

            # Check if we need to make that directory
            if not os.path.exists(directory):
                os.mkdir(directory)

            # Check if we have preexisting lightcurves
            if not os.path.exists('/'.join([directory, 'lightcurves'])):
                os.mkdir('/'.join([directory, 'lightcurves']))

        elif command == 'observatory':
            # Changes the observing location
            obsname = args[0]
            self.params['obsname'] = obsname
            printer("Observing location: {}".format(obsname))
        
        elif command == 'coords':
            # Changes the coordinates of the object you're about to talk about.
            if len(args) < 2:
                printer("I didn't get the right RA and Dec format! -> 'RA Dec'")
                printer("Please use:\n  RA - HH:MM:SS.SS\n  DEC - DD:MM:SS.SS\n")
                pass
            else:
                coords = '{} {}'.format(args[0], args[1])
                coords.replace(':', ' ')
                self.params['coords'] = coords
                printer("Using the following star coordinates:\n  RA:  {}\n  Dec: {}".format(args[0], args[1]))
        
        elif command == 'extinction':
            ext = float(args[0])
            self.params['ext'] = ext
            printer("Extinction coefficient: {}".format(ext))
        
        elif command == 'writeparams':
            if args == None:
                paramname = 'reduction_params.txt'
            else:
                paramname = args[0]
            with open(paramname, 'w') as f:
                for i, item in enumerate(self.params):
                    f.write("{} {}\n".format(item, self.params[item]))
            printer("Wrote parameters to 'reduction_params.txt'!")


        # SDSS field observations calibration
        elif command == 'sdss':
            # Toggle SDSS field
            if args == []:
                printer("Missing argument! Usage:\nSDSS [y/yes/True/n/no/False]")
            SDSS = args[0] in ['y', '1', 'yes', 'true']
            self.params['SDSS'] = SDSS
            printer("Are we in the SDSS field? [{}]".format(SDSS))
        
        elif command == 'stdcoords':
            # Changes the coordinates of the object you're about to talk about.
            if len(args) < 2:
                printer("I didn't seem to get the right RA and Dec format! -> 'RA Dec'")
                printer("Please use:\n  RA - HH:MM:SS.SS\n  DEC - DD:MM:SS.SS\n")
                pass
            else:
                coords = '{} {}'.format(args[0], args[1])
                coords.replace(':', ' ')
                self.params['stdcoords'] = coords
                printer("Using the following standard star coordinates:\n  RA:  {}\n  Dec: {}".format(args[0], args[1]))

        elif command == 'stdlogfile':
            if args == []:
                printer("Didn't get a file! Usage:\nstdLogfile [file]")
            stdLogfile = args[0]
            self.params['stdLogfile'] = stdLogfile
            printer("Using the standard star in this log file,\n  {}".format(stdLogfile))
        
        elif command == 'comparisonlogfiles':
            # Read in log filenames, terminated by an empty line, i.e. in the format:
            # logfiles
            # file1
            # file2
            # file3 
            #
            # <continue>
            fnames = []
            line = self.inFile.readline().strip()
            while line!='':
                if '#' in line: 
                    while line[0] == '#':
                        line = self.inFile.readline().strip()
                    if '#' in line:
                        line = line[:line.index('#')].strip()
                fnames.append(line)
                line = self.inFile.readline().strip()

            printer("Using the following logfiles for calculating the SDSS correction on each eclipse:")
            for fname in fnames:
                printer("- {}".format(fname))
            if fnames == []:
                printer("Didn't see any files! Files are read in line by line after the comparisonLogFiles command, and are terminated by a blank line.")
            printer("")

            self.params['comparisonfnames'] = fnames

        elif command == 'stdmags':
            # Must be in the format <r' g' u'>
            if len(args) < 3:
                printer("I didn't get enough magnitudes for the standard star!")
            else:
                mags = [float(m) for m in args[:3]]
                self.params['mags'] = mags

                printer("The standard star has the following apparent magnitudes:")
                printer("  r': {:.3f}\n  g': {:.3f}\n  u': {:.3f}".format(
                    mags[0], mags[1], mags[2]
                ))
        


        # getKappa stuff <-- Depricated!
        elif command == 'getkappa':
            self.getKappa()
        
        elif command == 'kappa_corr':
            kappas = [np.nan, np.nan, np.nan, np.nan]
            kappas[1:] = [float(x) for x in args[:3]]
            self.params['kappas']

            printer("Using the Kappa Corrections:")
            printer("  r': {}\n  g': {}\n  u': {}".format(
                kappas[1], kappas[2], kappas[3]
            ))



        # Eclipse times, and ephemeris stuff
        elif command == 'geteclipsetimes':
            self.getEclipseTimes()
        
        elif command == 'period':
            period = float(args[0])
            self.params['period'] = period
            printer("Using the period: {}".format(period))

        elif command == 't0':
            T0 = float(args[0])
            self.params['T0'] = T0
            printer("Using the T0: {}".format(T0))

        elif command == 'fitephemeris':
            self.fitEphem()


        # combineData
        elif command == 'combinedata':
            self.combineData()
        elif command == 'oname':
            if args == []:
                printer("Warning! Didn't get a filename!\nUsage: oname [file]")
            oname = args[0]
            self.params['oname'] = oname
            printer("Using the following filename: {}".format(oname))
        elif command == 'binsize':
            binsize = int(args[0])
            self.params['binsize'] = binsize
            printer("Binning data by {}".format(binsize))
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
                if '#' in line: 
                    while line[0] == '#':
                        line = self.inFile.readline().strip()
                    if '#' in line:
                        line = line[:line.index('#')].strip()
                fnames.append(line)
                line = self.inFile.readline().strip()

            self.params['fnames'] = fnames

            printer("Using the following logfiles:")
            for fname in fnames:
                printer("- {}".format(fname))
            printer("")

        # plotAll 
        elif command == 'overplot':
            if args != []:
                oname = args[0]
            else:
                oname = ''
            plot_all(self.written_files, oname)


        # Unknown command handler
        elif command not in ['', None]: # I don't want to be told about every blank line...
            printer("Unknown command!")
            printer("Command: {}".format(command))
            printer("   Args: {}".format(args))


if __name__ == "__main__":
    inf = None
    prom = []
    f = sys.argv

    if len(f) == 1:
        interp = Interpreter(prompt=['help'])

    if os.path.isfile(f[1]):
        infile = f[1]
        interp = Interpreter(inFile=infile)
    else:
        interp = Interpreter(prompt=True)
