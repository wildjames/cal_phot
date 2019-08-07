#!/usr/bin/env python3

import numpy as np
import sys
import os

from getEclipseTimes import getEclipseTimes
from combineData import combineData
from fitEphem import fitEphem
from plotAll import plot_all

from logger import printer, header

class Interpreter:
    def __init__(self, inFile=None, prompt=False):
        # Resize the terminal
        print(r"\x1b[8;{110};{80}t")
        print(r"\033[2J")

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

        self.instruments = ['hcam', 'ucam', 'uspec']

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
        -- ...
        -- Reduced_Data
        --- Night 1 reduction
        --- Night 2 reduction
        --- ...

        If you want to ensure the script works properly, it's best to follow this convention.

        The following are all case-insensitive. A hash causes the interpreter to ignore the rest of the line.

        - *CombineData*:
                Triggers the actual flux calibration script, using the supplied parameters.
        - *ComparisonLogFiles*:
                Starts reading the input file for the comparison star reductions, reduced with the
                same settings as the standard.
                Reads each following line for a filename, until it finds a blank line.
        - *Coords* [RA, str] [Dec, str]:
                RA and Dec of the target star
        - *Directory* [str]:
                Change the working location of the script. If it doesn't exist, create it.
        - *Extinction* [CCD1, float] [CCD2, float] [CCD3, float] ...:
                Extinction coefficients, in order of CCD. i.e. for ULTRACAM, r', g', 'u.
        - *fitEphemeris*:
                Take the ephemeris data contained in eclipse_times.txt (if this doesn't exist,
                creates it), and fits period and T0 to it starting with the previously found
                values.
        - *GetEclipseTimes*:
                Search the directory supplied by the *directory* command for .log files, and uses
                them to search for eclipse times. Saves to file
        - *Help*:
                Print the help string
        - *LogFiles*:
                Read in filenames for the target reduction. One file per line, terminated by an
                empty line.
        - *Observatory* [str]:
                Change the observation location. Must be interpreted by Astropy
                ([list here](http://docs.astropy.org/en/stable/api/astropy.coordinates.EarthLocation.html#astropy.coordinates.EarthLocation.of_site))
        - *oname* [str]:
                User-supplied template for writing the lightcurve files.
        - *overplot* [filename, str] [band, str]:
                plots the lightcurves from the given band (e.g. u, g, r) and saves to the filename
                given
        - *Period* [P, float]:
                The previously known period of the system
        - *SDSS* [bool]:
                Tell the script whether we're in the SDSS field or not. If we are, do a lookup for
                comparison magnitudes. If not, use a standard star observation.
        - *StdCoords* [RA, str] [Dec, str]:
                RA and Dec of the standard observation, if we're not in the SDSS.
        - *StdLogfile* [file, str]:
                Tells the script where to look for the standard star reduction
        - *StdMags* [CCD1, float] [CCD2, float] [CCD3, float]...:
                The apparent magnitude of the standard star, in CCD order.
        - *Stop/Exit/Quit*:
                Halt the reduction
        - *T0* [T0, float]:
                The previously known T0 of the system
        - *Writeparams* [filename, str]:
                Write out all the parameters as they stand at that point, to the file given
        ''')

    def getEclipseTimes(self):
        fnames    = self.params['fnames']
        coords    = self.params['coords']
        obsname   = self.params['obsname']
        directory = self.params['directory']

        printer("Getting eclipse times from data...")

        getEclipseTimes(fnames, coords, obsname, myLoc=directory)


    def fitEphem(self):
        directory = self.params['directory']
        T0        = self.params['T0']
        period    = self.params['period']

        printer("Fitting eclipse times to refine period and T0 parameters")

        T0, period = fitEphem(directory, T0, period)
        self.params['T0'] = T0
        self.params['period'] = period

    def combineData(self):
        oname     = self.params['oname']
        coords    = self.params['coords']
        obsname   = self.params['obsname']
        inst      = self.params['inst']
        T0        = self.params['T0']
        period    = self.params['period']
        myLoc     = self.params['directory']
        fnames    = self.params['fnames']
        SDSS      = self.params['SDSS']
        ext       = self.params['ext']

        printer("Combining, calibrating, and plotting data...")
        if SDSS:
            written_files = combineData(oname, coords, obsname, T0, period, SDSS=True, inst=inst,
            myLoc=myLoc, fnames=fnames, ext=ext)
        else:
            # Retrieve the SDSS-matching reductions for each night
            comparisons = self.params['comparisonfnames']
            stdLogfile  = self.params['stdLogfile']
            stdCoords   = self.params['stdcoords']
            stdMags     = self.params['mags']

            # ref_kappa = self.params['kappas']
            written_files = combineData(oname, coords, obsname, T0, period, SDSS=False,
                myLoc=myLoc, fnames=fnames, comp_fnames=comparisons, inst=inst,
                std_fname=stdLogfile, std_coords=stdCoords, std_mags=stdMags, ext=ext)

        self.written_files += written_files

        self.written_files = sorted(self.written_files)

        printer("So far, I've written the following files:")
        for f in self.written_files:
            printer("-> {}".format(f))

    def parse(self, line):
        line = line.strip()

        # Clean up input for parsing
        if '#' in line:
            cut = line.index('#')
            line = line[:cut]

        if line == '':
            command = None
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
        if command == 'help':
            self.help()
            exit()
        elif command in ['stop', 'exit', 'quit']:
            printer("Bye!")
            print("\n\n")
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
            lc_dir = os.path.join(directory, 'MCMC_LIGHTCURVES')
            if not os.path.exists(lc_dir):
                os.mkdir(lc_dir)

        elif command == 'observatory':
            # Changes the observing location
            obsname = ' '.join(args)
            self.params['obsname'] = obsname
            printer("Observing location: '{}'".format(obsname))

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
            ext = [float(i) for i in args]
            self.params['ext'] = ext
            printer("Extinction coefficients (in order of CCD): {}".format(ext))

        elif command == 'writeparams':
            if args == None:
                paramname = 'reduction_params.txt'
            else:
                paramname = args[0]

            with open(paramname, 'w') as f:
                for item in self.params:
                    f.write("{} {}\n".format(item, self.params[item]))
            printer("Wrote parameters to 'reduction_params.txt'!")

        elif command == 'inst':
            if args == None:
                raise Exception("You need to define an instrument!")
            elif args[0] in self.instruments:
                self.params['inst'] = args[0]
                print("Observations were taken with {}".format(self.params['inst']))

        # SDSS field observations calibration
        elif command == 'sdss':
            # Toggle SDSS field
            if args == []:
                printer("Missing argument! Usage:\nSDSS [y/yes/True/1/0/n/no/False]")
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
                while '#' in line:
                    # If the hash is not in the first space, cut the line to size
                    if '#' in line[1:]:
                        line = line[:line.index('#')].strip()
                    # If the hash is in the first space, read the next line
                    else:
                        line = self.inFile.readline().strip()
                if line != '':
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
                while '#' in line:
                    # If the hash is not in the first space, cut the line to size
                    if '#' in line[1:]:
                        line = line[:line.index('#')].strip()
                    # If the hash is in the first space, read the next line
                    else:
                        line = self.inFile.readline().strip()
                if line != '':
                    fnames.append(line)
                line = self.inFile.readline().strip()

            self.params['fnames'] = fnames

            printer("Using the following logfiles:")
            for fname in fnames:
                printer("- {}".format(fname))
            printer("")

        # plotAll
        elif command == 'plot' or command == 'overplot':
            if args != []:
                oname = args[0]
                band  = args[1]
            else:
                oname = ''
            myLoc = self.params['directory']
            plot_all(self.written_files, oname, band, myLoc)


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
        interp = Interpreter(prompt=True)

    if os.path.isfile(f[1]):
        infile = f[1]
        interp = Interpreter(inFile=infile)
