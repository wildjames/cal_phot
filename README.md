# cal_phot
A set of scripts that will calibrate periodic \*CAM lightcurves from counts to flux above atmosphere.

The scripts available are fairly basic calibration stuff, but this takes a lot of the legwork out of reducing periodic data. Importantly, it makes calls to SDSS, either for the comparisons in the target frame, or a separate reference star taken on the same night/run. And don't worry, these are all corrected for airmass, light travel time, etc! The overall output are a few files:
- A lightcurve for the star in aperture 1 of each supplied `.log` file
- Two figures for each lightcurve, one showing the target's flux over the observation, and another plotting the comparisons against each other
- An mcmc chain that optimises the ephemeris data

The figure containing the compared comparison stars is there to help diagnose a bad reference star. These should all ideally be flat lines at unity, but life is rarely that perfect.

The main advantage of using this script, and in writing this interpreter at all, is that it creates a consistent record of how a data set was reduced for future reference. This will reduce 'black box' complaints later, when it comes to writing up what we find in the lightcurve data.

## bin_data.py
Additionally, this git contains my binning script. This script makes an evenly spaced 1D grid of phase times, and combines the relevant data from given files into those bins. It respects errors, weights means by them, and uses the mean phase of the binned data in each bin to preserve as much information a possible.

# Using the Pipeline
This was written as a side-project, and hasn't been streamlined very much. However, there is a recommended workflow. The process is slightly different for stars in or out of the SDSS field.

## Stars in the SDSS
### 1. Reduce your raw data with the HIPERCAM pipeline
The HiPERCAM pipeline is fairly easy to use. There's great doumentation to be found [here](http://deneb.astro.warwick.ac.uk/phsaap/hipercam/docs/html/commands.html), but I'll summarise the rough steps;

- Make a mean frame, used to define apertures
  - `averun`
- Set the apertures that you want to extract
  - `setaper`
- Generate the reduction file
  - `genred`
- Tweak the reduction file and optimise the settings
  - (black magic)
- Reduce the data!
  - `reduce`

This produces a `.log` file, containing the electron counts per aperture, per frame. However, this still needs to be calibrated, and DOESN'T contain any pointing info. So, we need to do a lookup for some actual apparent magnitudes.


### 2. Make auxilliary files

I find it easiest to start thinking about claibration here, and later on `cal_phot` is going to want to know your comparison star magnitude. Create a new text file with the same name as your log file, and the extension `.coords`, i.e. I reduce a file and call it `star_data.log` with 2 comparison apertures, then create another file called `star_data.coords`. This file is going to contain the RA and Dec of our comparison stars, and the filter that we want the observations in. This is the general format:
```
[CCD1 filter] [CCD2 filter] [CCD3 filter] ...

[CCD1 AP1 RA] [CCD1 AP1 DEC]

[CCD2 AP1 RA] [CCD2 AP1 DEC]
[CCD2 AP2 RA] [CCD2 AP2 DEC]

[CCD3 AP3 RA] [CCD3 AP3 DEC]

...
```
Note that the whitespace is important, the apertures of the different CCDs are separated by a blank line. I find it easiest to have an `Aladin` window open on the side, where right-clicking on a star and clicking `Copy recticle position to clipboard` gives you the RA and Dec.

You should now have two files for each reduction, a `.log` and a `.coords`

### 3. Write an input file

This is the easy part. Create a raw text file, and enter the commands you want to run. They are exectured in order, so multiple instruments, or datasets, or even objects can be calibrated in one run. Though, for ease of understanding what's gone on, it's better to fragment your files. For a star in SDSS, the input is likely to look something like this:

```
# Where are we observing from?
observatory 18:35:26, +98:29:12
inst uspec

### Prior knowledge of the target system
coords      07:48:59.55 +31:25:12.6    # RA, Dec of target
period      0.0583110795
T0          57808.63029759

### Can we do an SDSS lookup? ###
SDSS 1

### What files contain the targets? ###

## This is a list of the best eclipses we have, for the first round of fitting. ##
logfiles    # List of logfiles to calibrate
2017-01-22_KG5.log
2017-02-15_KG5.log
2017-12-12_KG5.log

## Refine our ephemeris for the system ###
getEclipseTimes     # Get eclipse timing from data files
fitEphemeris        # Fit the ecplise times for better ephemeris
   
### Process the files ###
oname       SDSSJ0748_0
combineData

### Now I want to plot my results.
# overplot [filename] [pattern]
plot    SDSSJ0748_0_KG5   _KG5
```

That looks like a lot, but fairly easy to break down. First, we tell the script where the observations it's about to process were taken from, in `lat, lon`. This allows the code to correct our observations to heliocentric time, and remove slight timing differences due to the position of the earth changing throughout the year. The code also needs to know the RA and Dec of the target, for the same reason. In a non-SDSS calibration, this will also be used to calculate and subtract airmass effects.

We tell the script what log files contain the data we're interested in. The list is terminated by a blank line.

Then, we define our existing ephemeris information *in days* before running a script to extract the eclipse times from the `.log` files we've defined as containing eclipses. Once that's done, we tell the interpreter to run the ephemeris refining tool. This fits the eclipse time to all the lightcurves it finds in the working directory, but asks for confirmation before each one.

`oname` defines an output filename template, and combinedata pulls the trigger on the final processing. The processing is interactive and should be carefully monitored, but should you doubt yourself later, a log of the reduction will be generated in `Calibration.log`. Finally, I tell `cal_phot` to `plot` the data in a `.pdf`.

### 3. Run it!
`python3 interpreter.py commandFile.dat`

# TODO:
- Move the input file over to YAML, and make the input method less horiffic to use.


# COMMANDS
The following are all case-insensitive.

- A hash causes the interpreter to ignore the rest of the line

- *CombineData*:
  - Triggers the actual flux calibration script, using the supplied parameters.
- *ComparisonLogFiles*:
  - Starts reading the input file for the comparison star reductions, reduced with the same settings as the standard.
  - Reads each following line for a filename, until it finds a blank line.
- *Coords* \[RA, str] \[Dec, str]: 
  - RA and Dec of the target star
- *Directory* [str]: 
  - Change the working location of the script. If it doesn't exist, create it.
- *Extinction* \[CCD1, float] \[CCD2, float] \[CCD3, float] ...: 
  - Extinction coefficients, in order of CCD. i.e. for ULTRACAM, r', g', 'u.
- *fitEphemeris*:
  - Take the ephemeris data contained in eclipse_times.txt (if this doesn't exist, creates it), and fits period and T0 to it starting with the previously found values.
- *GetEclipseTimes*:
  - Search the directory supplied by the *directory* command for .log files, and uses them to search for eclipse times. Saves to file
- *Help*: 
  - Print the help string
- *inst* \[str]
  - Sets the instrument used for the observations. One of \[uspec, ucam, hcam]
- *LogFiles*:
  - Read in filenames for the target reduction. One file per line, terminated by an empty line.
- *Observatory* [str]: 
  - Change the observation location. Must be interpreted by Astropy ([list here](http://docs.astropy.org/en/stable/api/astropy.coordinates.EarthLocation.html#astropy.coordinates.EarthLocation.of_site))
- *oname* \[str]:
  - User-supplied template for writing the lightcurve files.
- *overplot* \[filename, str] \[band, str]:
  - plots the lightcurves from the given band (e.g. u, g, r) and saves to the filename given
- *Period* \[P, float]:
  - The previously known period of the system
- *SDSS* \[bool]:
  - Tell the script whether we're in the SDSS field or not. If we are, do a lookup for comparison magnitudes. If not, use a standard star observation.
- *StdCoords* \[RA, str] \[Dec, str]:
  - RA and Dec of the standard observation, if we're not in the SDSS.
- *StdLogfile* \[file, str]:
  - Tells the script where to look for the standard star reduction
- *StdMags* \[CCD1, float] \[CCD2, float] \[CCD3, float]...:
  - The apparent magnitude of the standard star, in CCD order.
- *Stop/Exit/Quit*: 
  - Halt the reduction
- *T0* \[T0, float]:
  - The previously known T0 of the system
- *Writeparams* \[filename, str]: 
  - Write out all the parameters as they stand at that point, to the file given


## Usage
- getEclipseTimes:
  - coords: Needed to correct the observation time to Barycentric MJD
  - obsname: Observing location, for the same
  - directory: The directory to search for .log files to analyse. Also saves eclipse times to this folder.
- fitEphem:
  - directory: The directory to search for prior eclipse times. Also saves the resulting corner plot from the fit.
  - T0: Initial value
  - Period: Initial value
- combineData:
  - oname: Naming template for created files.
  - coords: RA and Dec of the target star.
  - Observatory: Observatory location
  - T0: Ephemeris data
  - Period: Ephemeris data
  - directory: Location to put created files. If it's 'Reduced_Data/\*', put it there, if not, create a subfolder called 'Reduced_Data' and puts stuff in that.
  - logFileNames: Location of the target systems' .log files. 
  - SDSS: Are we in the SDSS field?
  - extinction: Extinction coeffiecients. 
