# GalSpec package
This package is intended to quickly and easily generate a galaxy spectrum with a blackbody continuum emission and emission lines. The lines in this package are:
* CO lines 
* SIII
* SiII
* OIII, OI 
* NIII, NII 
* CII, CI

*Code by Tom Bakx, packaging and addition of molecular lines and rescaling by Stefanie Brackenhoff*

## Functionalities & Usage
* A spectrum can be generated using the ```spectrum()``` function. 

It takes the following inputs:
* **luminosity** in units of log(L_fir [L_sol])
* **redshift** z
* **fLow, fHigh** minimum and maximum frequency that will be in the spectrum in units of GHz
* **numFreqBins** amount of linearly spaced frequency bins at which the spectrum should be evaluated
* **linewidth** width of spectral lines in units of km/s
* **COlines** 'Kamenetzky' [1] or 'Rosenberg' [2] determines the amplitude of the CO lines, default is Kamenetzky. 
* **lines** 'Bonato' [3] or 'Spinoglio' [4] determines the amplitude of the remaining spectral lines, default is Bonato
* **mollines** 'True' or 'False'. Toggles whether molecularlines are shown. Values estimated using [5]. Defaults as 'True'
* **variance** adds uncertainty to the amplitude of the atomic lines. Defaults to 0. Variation in molecular lines not implemented in the current version.
* **manualrescale** sets whether the lines should be default ratios by Kamenetzky/Rosenberg and Bonato/Spinoglio, or set by user. Options: 
	- **'False'** default line amplitudes used
    - **'Absolute'** numpy array of 37 line amplitudes in Jy additive to blackbody emission can be set in rescaleArray. Additional entries ignored if mollines is set to 'False'.
    - **'Relative'** numpy array of 37 scalars can be set in rescaleArray. The default ratios from Kamenetzky/Rosenberg are multiplied by these scalars prior to addition to the blackbody spectrum. Additional entries ignored if mollines is set to 'False'.
* **rescaleArray** rescales the emission lines according to setting in 'manualrescale'. Order of the lines can be found using ```linenames()```

And creates as output:
* **freqArray** array frequencies in units of GHz
* **spectrum** array of the flux densities in the spectrum in units of Jy

[1]: http://dx.doi.org/10.3847/0004-637X/829/2/93
[2]: https://ui.adsabs.harvard.edu/link_gateway/2015ApJ...801...72R/doi:10.1088/0004-637X/801/2/72
[3]: https://ui.adsabs.harvard.edu/link_gateway/2014MNRAS.438.2547B/doi:10.1093/mnras/stt2375
[4]: https://ui.adsabs.harvard.edu/link_gateway/2012ApJ...745..171S/doi:10.1088/0004-637X/745/2/171
[5]: https://arxiv.org/pdf/1106.5054.pdf



* The spectrum can quickly be plotted using the ```plotspectrum()``` function

This function takes the outputs of ```spectrum()``` as an input and creates a plot with axis labels



* The names of the spectral lines and their order for the ```rescaleArray``` in ```spectrum()``` are outputted in a numpy array. Whether the names of molecular lines are shown can be toggled by setting the keyword mollines to 'False'.

## Examples
* Simple example
```
import galspec

luminosity = 13.7
z = 4.43
fLow = 100 #GHz
fHigh = 900 #GHz
numFreqBins = 1500
linewidth = 600
gal_freq, gal_flux = galspec.spectrum(luminosity, z, fLow, fHigh, numFreqBins, linewidth, mollines = 'False')

galspec.plotspectrum(gal_freq, gal_flux)
```

![Example](/example_spectrum.png)

* Using the ratios by Rosenberg and Spinoglio
```
import galspec

luminosity = 13.7
z = 4.43
fLow = 100 #GHz
fHigh = 900 #GHz
numFreqBins = 1500
linewidth = 600
gal_freq, gal_flux = galspec.spectrum(luminosity, z, fLow, fHigh, numFreqBins, linewidth, 'Rosenberg', 'Spinoglio', mollines = 'False')

galspec.plotspectrum(gal_freq, gal_flux)
```
![Example2](/example_spectrum_RosSpin.png)

* Using manual ratios to only show CO lines
```
import galspec
import numpy as np

luminosity = 13.7
z = 4.43
fLow = 100 #GHz
fHigh = 900 #GHz
numFreqBins = 1500
linewidth = 600

names = galspec.linenames()
lines = np.zeros(len(names))
for i in range(len(names)):
    if names[i].startswith('CO'): lines[i]=1

gal_freq, gal_flux = galspec.spectrum(luminosity, z, fLow, fHigh, numFreqBins, linewidth, manualrescale = 'Relative', rescaleArray = lines, mollines = 'False')

galspec.plotspectrum(gal_freq, gal_flux)
```
![Example3](/example_spectrum_CO_only.png)

* Including molecular lines
```
import galspec
import numpy as np

luminosity = 13.7
z = 4.43
fLow = 100 #GHz
fHigh = 900 #GHz
numFreqBins = 1500
linewidth = 600


gal_freq, gal_flux = galspec.spectrum(luminosity, z, fLow, fHigh, numFreqBins, linewidth)

galspec.plotspectrum(gal_freq, gal_flux)
```
![Example4](/example_spectrum_with_mols.png)

## Installation
```pip install galspec```

## Required packages
* ```Numpy```
* ```astropy```
* ```matplotlib```
* ```scipy```