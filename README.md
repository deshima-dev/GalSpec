# GalSpec package
This package is intended to quickly and easily generate a galaxy spectrum with a blackbody continuum emission and CO and CII emission lines.
Code by Tom Bakx, packaging by Stefanie Brackenhoff

## Functionalities & Usage
* A spectrum can be generated using the ```spectrum()``` function. 

It takes the following inputs:
**luminosity** in units of log(L_fir [L_sol])
**redshift** z
**fLow, fHigh** minimum and maximum frequency that will be in the spectrum in units of GHz
**numFreqBins** amount of linearly spaced frequency bins at which the spectrum should be evaluated
**linewidth** width of spectral lines in units of km/s

And creates as output:
**freqArray** array frequencies in units of GHz
**spectrum** array of the flux densities in the spectrum in units of Jy

* The spectrum can quickly be plotted using the ```plotspectrum()``` function
This function takes the outputs of ```spectrum()``` as an input and creates a plot with axis labels

## Example
```
import galspec

luminosity = 13.7
z = 4.43
fLow = 332 #GHz
fHigh = 377 #GHz
numFreqBins = 1500
linewidth = 600
gal_freq, gal_flux = galspec.spectrum(luminosity, z, fLow, fHigh, numFreqBins, linewidth)

galspec.plotspectrum(gal_freq, gal_flux)
```

![Example](/example_spectrum.png)


## Installation
```pip install galspec```

## Required packages
* ```Numpy```
* ```astropy```
* ```matplotlib```