# *CAM flux calibration pipeline

A set of scripts that will calibrate periodic \*CAM lightcurves from counts to flux above atmosphere.

The scripts available are fairly basic calibration stuff, but this takes a lot of the legwork out of reducing periodic data. Importantly, it makes calls to SDSS, either for the comparisons in the target frame, or a separate reference star taken on the same night/run. And don't worry, these are all corrected for airmass, light travel time, etc! The useful outputs are a few files:

-   A lightcurve for the star in aperture 1 of each supplied `.log` file
-   Two figures for each lightcurve, one showing the target's flux over the observation, and another plotting the comparisons against each other
-   An mcmc chain that optimises the ephemeris data

The figure containing the compared comparison stars is there to help diagnose a bad comparison. These should all ideally be flat lines at unity, but life is often not that perfect.

The main advantage of using this script, and in writing this at all, is that it creates a consistent record of how a data set was reduced for future reference. This will reduce 'black box' complaints later or, more relevant to me, dodge the problem of forgetting if I did it 100% right, when it comes to writing up what we find in the lightcurve data.

A fairly verbose log is produces by each script in the file, `Calibration.txt`. Note that if an older one exists where it wants to put a new one, it'll overwrite the old log. Normally you'll want this - so that your logs always correspond to the data you currently have, but please beware!

# Installing the package

This package is pip installable:

```
pip3 install calphot
```

but if you want to do it more manually, and fiddle with the source:
```
git clone https://github.com/wildjames343/cal_phot
cd cal_phot
pip3 install --user .
```

This should also create some commands for you to use from the command line;
  - `cal_phot`
  - `calc_extinction`
  - `comparison_mags`
  - `bin_data`

All have `-h` doc strings, and `cal_phot` has two [example](example_files/cal_commands_SDSS.yaml) [files](example_files/cal_commands_noSDSS.yaml)

# Scripts you'll want to use

## `cal_phot`

This is the main calibration script. See below for details, but in short it needs a configuration yaml file.

## `calc_extinction`

This takes a long observations and, assuming the apertures you use to reduce it are around constant-brightness stars, calculates the extinction in each frame. This requires a few things to be thought about to work properly:

  - You need to use an observation on a **clear, photometric** night. Otherwise, this will be wrong! Clouds would contribute to the extinction and give you a wrong value.
  - The larger the airmass range, the better. The logs on `deneb` list the airmass range for an observation, so use that to inform what run you use.
  - Variable targets will result in junk. There's a lazy way and a smart way to deal with this - use many apertures to make sure that at least some are constant sources, or actually check each target's RA and Dec in catalogues for variability.

## `comparison_mags`

We need to know the apparent magnitudes of the comparison stars so that we can convert ADU counts to a flux in mJy. This is easy when we're in the SDSS field, since the code can just look that up, but if we're elsewhere in the sky we have to jump through some hoops. This script will calculate the average apparent magnitude of all apertures in a logfile, except the first one (assumed to be the target). The result of this can be put into the `cal_phot` config file.

This should be done carfully - there's no point in flux calibrating the target against something that's not true! Take an observation from a clear, photometric night, and use this to calculate the magnitudes of ALL the comparison stars you use - if some are calculated on a night with even a little cloud, that magnitude will be WRONG!

## `bin_data`

Note to self: **Always be wary about *how* and *why* you're binning some data together!**

This git contains my binning script, which takes a few phase-folded lightcurves and averages them together.

This script makes an evenly spaced 1D grid of phase times, and combines the relevant data from given files into those bins. It respects their errors *and* weights data by them, and uses the mean phase of the binned data in each bin to preserve as much information as possible.

# Using my Pipeline

This was written as a side-project, and hasn't been streamlined very much. However, I do have a recommended workflow. The process is slightly different for stars in or out of the SDSS field.

# Stars in the SDSS

## 1. Reduce your raw data with the HIPERCAM pipeline

The HiPERCAM pipeline is fairly easy to use. There's great doumentation to be found [here](http://deneb.astro.warwick.ac.uk/phsaap/hipercam/docs/html/commands.html), but I'll summarise the rough steps;

-   Make a flat and bias frame from the relevant runs
    -   `makeflat`, `makebias`
-   Make a mean frame, used to define apertures
    -   `averun`
-   Set the apertures that you want to extract
    -   `setaper`
-   Generate the reduction file
    -   `genred`
-   Tweak the reduction file and optimise the settings
    -   (((black magic)))
-   Reduce the data!
    -   `reduce`

This produces a `.log` file, containing the electron counts per aperture, per frame. However, this still needs to be calibrated, and DOESN'T contain any pointing info! We need to do a lookup for some actual apparent magnitudes.

## 2. Make `.coords` files

Create a new text file with the same name as your log file, and the extension `.coords`, i.e. I reduce a file and call it `star_data.log` with 2 comparison apertures, then create another file called `star_data.coords`. This file is going to contain the RA and Dec of our *comparison* stars (**there is no entry for the target!!**), and the filter that we have the observations in. This is the general format:

```
    [CCD1 filter] [CCD2 filter] [CCD3 filter]

    [CCD1 AP1 RA] [CCD1 AP1 DEC]
    [CCD1 AP2 RA] [CCD1 AP2 DEC]

    [CCD2 AP1 RA] [CCD2 AP1 DEC]
    [CCD2 AP2 RA] [CCD2 AP2 DEC]

    [CCD3 AP2 RA] [CCD3 AP2 DEC]
    [CCD3 AP3 RA] [CCD3 AP3 DEC]
```

Note that the whitespace is important, the apertures of the different CCDs are separated by a blank line. I find it easiest to have an `Aladin` window open on the side, where right-clicking on a star and clicking `Copy recticle position to clipboard` gives you the RA and Dec in the right format. `cal_phot` will use these to get those stars' magnitudes from the SDSS.

**You should now have two files for *each* reduction**, a `.log` and a `.coords`. If a `coords` file is missing, `cal_phot` will let you know about it.

## 3. Calculate the extinction coefficients of the night

The SDSS standard and the target frame are, by definition, oberved at different locations in the sky. If they weren't, we wouldn't need the standard! So, we need to correct for atmospheric extinction. The script, `calc_extinction`, helps with this. Save the output of this for the input file. This isn't always 100% necessary, as the observatories each have fairly consistent values that they publish.

## 4. Write an input file

This is the easy part. The config is a YAML file, and I have an example of what it should look like in this git, [here!](cal_commands_SDSS.yaml)

## 5. Run it!

`cal_phot cal_commands.yaml`

# Non-SDSS Systems

These are slightly trickier, since I can't do an automatic lookup for `*CAM`-like filters. You'll have to help me out with a bit of extra legwork by also reducing a standard star observation on (or near) the night that the target was observed. The script is then going to take that observation, knowing the magnitude of the star, and work out the electron-count-to-flux ratio which can then be applied to the target reductions.

## 1. Reduce the target observations

This is the same as for an SDSS-field system. Just do as you would normally.

## 2. Reduce the standard star observation

The observer has likely observed a standard on the night (check the logs). Reduce this system with settings that make sense for it, but REMEMBER THOSE SETTINGS! In order to be consistent between the std. and the comparisons that we care about later, they must be reduced with all the same settings. Also use a large, fixed, aperture size to be sure that you're not getting *any* light leaking out of the edge of the target aperture.

Then, go over to the target obeservation, and use the **exact same settings** to reduce the target frames. *Use the same `.ape` file as you used in step 1*, but make sure you use the same settings that you have in the standard star `.red` file! This will ensure consistency. I tend to call this standard-like reduction of the target frame `<system_name>_standard.log`.

## 3. Compute the comparison star magnitudes

Use the script `comparison_mags.py` to get the apparent magnitudes of all your comparison stars. This script will output the magnitudes in a conviniently copy-pasteable chunk for your configuration file.

## 4. Construct the configuration file

There are a few extra bits of info needed in the configuration now. I usually just take the [example](cal_commands.yaml), copy it into my working directory, and tweak the stuff.


## 5. Run the interpreter

That's it! The software will walk you through the rest.

# Problems?

Email me, so I can find out what I've explained badly and improve this walkthrough.

* * *
