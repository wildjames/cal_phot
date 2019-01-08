##TODO: I need to rework the handling of myLoc, and the directory command

# cal_phot
A set of scripts that will calibrate *CAM lightcurves from counts to flux above atmosphere.

These are currently fairly inflexible. Currently, this is a list of this packages' functions:
- Given a hipercam pipeline's .log file for a star, fold and plot it about its ehemeris.
- If the star is in the SDSS field, perform a lookup for given comparison stars and convert the target's photon counts into apparent magnitude above the atmosphere.
- Given the SDSS magnitude of a standard star, and its own separately observed .log file, convert non-SDSS field observations to apparent magnitudes above the atmosphere.
- Fit ecplise times of a set of lightcurves, and use this to calculate the ephemeris data for a system.

The advantage of using this script, and in writing this interpreter at all, is that it creates a consistent record of how a data set was reduced for future reference. This will reduce 'black box' complaints later, when it comes to writing up what we find in the lightcurve data.

# COMMANDS
The following are all case-insensitive.

- A hash causes the interpreter to ignore the rest of the line
- *Help*: 
  - Print the help string
- *Stop/Exit/Quit*: 
  - Halt the reduction
- *Directory* [str]: 
  - Change the working location of the script. If it doesn't exist, create it.
- *Observatory* [str]: 
  - Change the observation location. Must be interpreted by Astropy ([list here](http://docs.astropy.org/en/stable/api/astropy.coordinates.EarthLocation.html#astropy.coordinates.EarthLocation.of_site))
- *Coords* \[RA, str] \[Dec, str]: 
  - RA and Dec of the target star
- *Extinction* \[CCD1, float] \[CCD2, float] \[CCD3, float] ...: 
  - Extinction coefficients, in order of CCD. i.e. for ULTRACAM, r', g', 'u.
- *Writeparams* \[filename, str]: 
  - Write out all the parameters as they stand at that point, to the file given
- *Sdss* \[bool]:
  - Tell the script whether we're in the SDSS field or not. If we are, do a lookup for comparison magnitudes. If not, use a standard star observation.
- *StdCoords* \[RA, str] \[Dec, str]:
  - RA and Dec of the standard observation, if we're not in the SDSS.
- *StdLogfile* \[file, str]:
  - Tells the script where to look for the standard star reduction
- *ComparisonLogFiles*:
  - Starts reading the input file for the comparison star reductions, reduced with the same settings as the standard.
  - Reads each following line for a filename, until it finds a blank line.
- *StdMags* \[CCD1, float] \[CCD2, float] \[CCD3, float]...:
  - The apparent magnitude of the standard star, in CCD order.
- *GetEclipseTimes*:
  - Search the directory supplied by the *directory* command for .log files, and uses them to search for eclipse times. Saves to file
- *Period* \[P, float]:
  - The previously known period of the system
- *T0* \[T0, float]:
  - The previously known T0 of the system
- *fitEphemeris*:
  - Take the ephemeris data contained in eclipse_times.txt (if this doesn't exist, creates it), and fits period and T0 to it starting with the previously found values.
- *CombineData*:
  - Triggers the actual flux calibration script, using the supplied parameters.
- *oname* \[str]:
  - User-supplied template for writing the lightcurve files.
- *binsize* \[int]:
  - CURRENTLY REMOVED. DO MANUALLY AFTER THE FACT. 
- *LogFiles*:
  - Read in filenames for the target reduction. One file per line, terminated by an empty line.
- *overplot* \[filename, str] \[band, str]:
  - plots the lightcurves from the given band (e.g. u, g, r) and saves to the filename given
