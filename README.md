# cal_phot
A set of scripts that will calibrate *CAM lightcurves from counts to flux above atmosphere.

These are currently fairly inflexible. Currently, this is a list of this packages' functions:
- Given a hipercam pipeline's .log file for a star, fold and plot it about its ehemeris.
- If the star is in the SDSS field, perform a lookup for given comparison stars and convert the target's photon counts into apparent magnitude above the atmosphere.
- Given the SDSS magnitude of a standard star, and its own separately observed .log file, convert non-SDSS field observations to apparent magnitudes above the atmosphere.
- Fit ecplise times of a set of lightcurves, and use this to calculate the ephemeris data for a system.

The results of each are written to a file, except the magnitude corrections. [TODO: add this!]

The advantage of using this script, and in writing this interpreter at all, is that it creates a record of how a data set was reduced for future reference. This will reduce 'black box' complaints later, when it comes to writing up what we find in the lightcurve data.

This is a copy of the 'help' output from the script. This can also be given by running `interpreter.py help`

### COMMANDS

## Housekeeping commands
- `test <args>` -- Testing command, literally just prints the arguments.
- `help` -- Print this text
- `stop/exit/quit` -- Stops the processing
- `observatory <name>` -- Sets the observing location to the argument it's given. Run astropy.EarthLocation.get_site_names() to get a list of valid sites.
- `coords <RA> <Dec` -- Sets the sky location. This has to be changed if your standard is in a different place than your target!
- `directory <dir>` -- If we are using sub-folders, this tells the next script where to look for data
- `extinction <ext> -- Sets extinction coefficient, in mags/airmass
- `writeparams` -- Writes parameters to file. You probably want this every time!

## Kappa corrections commands
- `stdLogfile <file>` -- Tells the kappa correction script where to look for a .log file
- `stdMags <r' mag> <g' mag> <u' mag>` -- SDSS magnitudes for the standard star
- `kappa_corr <r' correction> <g' correction> <u' correction>` -- Manually define the correction factors

## Ephemeris commands
- `period <P>` -- Sets the period of the system for folding. If we are doing getEclipseTimes, this forms the initial guess
- `T0 <T0>` -- Sets the T0 for folding. 

## Data combination commands
- `oname <String>` -- Define the names prefix to use when writing out reduced data
- `binsize <int>` -- Define the binning after we've folded our data. <1> sets to unbinned.

## 'Trigger' commands
- `getKappa` -- Calls the script to calculate the photometric -> magnitude corrections.
    Requires:
     - An hcam logfile containing the photometry of the standard star.
     - Std. star coordinates
     - Std. star SDSS magnitudes
     - An extinction coefficient, in mags/airmass. This defaults to 0.161
     - The name of the observatory where the data were taken. Used to lookup in astropy
  
- `getEclipseTimes` -- Calls a script that searches the given directory for logfiles, and with the user's help determines the eclipse times of the lightcurves. These are then used to compute the best period for the data.
    Requires: 
     - Target coordiates
     - Initial period guess
     - [OPTIONAL] Previously found value of T0

- `combineData` -- Searches for logfiles, callibrates the photometry, folds and bins the data.
    Requires: 
     - Target coordinates
     - Kappa corrections
     - Observatory name
     - Ephemeris data (T0 and period)
     - Binsize
     - Directory
     - Name template to use while writing out to files
