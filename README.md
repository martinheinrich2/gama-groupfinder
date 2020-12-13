# gama-groupfinder

This is a number of scrips to find groups of galaxies from a volume-limited
sample of the GAMA-Survey. The volume limited sample was created from two GAMA-Survey tables.

The main data http://www.gama-survey.org/dr3/data/cat/EqInputCat/v46/TilingCat.fits
and distances from
http://www.gama-survey.org/dr3/data/cat/LocalFlowCorrection/v14/DistancesFrames.fits
were combined and a volume-limited sample was created. This was done with TOPCAT (http://www.star.bristol.ac.uk/~mbt/topcat/).

The approach is similar to a friends-of-friends algorithm,
here the steps are split in separate scrips.
The idea was to find groups for different
linking-lengths and to determine if the group distributions indicate a
special linking-length.

## Documentation

## Installation

## Basic features

1. Use Distance_calculator.py to calculate the line-of-sight and transverse
   distances between galaxies for a distinct cosmology.

2. Use galaxypair_finder.py to find galaxy pairs within a given range
   of line-of-sight (in km/s) and transverse (in Mpc) values.

3. Use galaxygroup_finder.py to find galaxy groups from galaxy pairs

4. Use fit_distribution_powerlaw.py to calculate powerlaw-fit and parameters

## Ideas

- It might be possible to create a distance matrix and use it in step 2. Depending on the
available memory of the computer, a distance matrix of 56119 galaxies is probably not easy to
handle.
