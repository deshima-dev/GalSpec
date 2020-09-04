##-----------------------------------------------
## Code written and edited by Tom Bakx
## Converted to a package by Stefanie Brackenhoff
## tjlcbakx@gmail.com
##-----------------------------------------------


##-----------------------------------------
## Header imports & colours
##-----------------------------------------

import numpy as np
from astropy.cosmology import Planck15 as cosmo
from matplotlib import pyplot as plt

orange = '#ff9500'#(1,0.584,0)
blue =  '#007aff'  #(0,.478,1) blue
green = '#4cd964'
red = '#ff3b30'
grey = '#8e8e93'   #(142./255,142./255,147./255)

##-----------------------------------------
## Necessary constants
##-----------------------------------------
T_cold     = 20.379
T_hot      = 43.7
Ratio       = 28.67
Beta        = 1.97

##-----------------------------------------
## Fit a bb-template for redshift and amplitude
##-----------------------------------------

def tomModel(v,Amp,z,T_cold,T_hot,Ratio,Beta):
    # Give the flux densities of the modeled SED at the requested frequencies v, at redshift z
    v_rest = (1+z)*v
    return Amp*((blackBody(v_rest,T_hot)*(v_rest**Beta)) + Ratio*(blackBody(v_rest,T_cold)*(v_rest**Beta)))

def blackBody(v,T):
    h = 6.626e-34
    c = 3.0e8
    k = 1.38e-23
    from numpy import exp
    return (2*h*(v*v*v))/(c*c*exp((h*v)/(k*T)) - 1)

def giveAmplitude(model,measurement,uncertainty,sumsquaresum=False):
    # This calculates the amplitude algebraically, to speed up computation time by
    # Decreasing the number of variables by the number of observations
    var1 = 0.
    var2 = 0.
    if sumsquaresum:
        for i in range(model.shape[0]):
            var1 += measurement[i]/uncertainty[i]
            var2 += model[i]/uncertainty[i]
    else:
        for i in range(model.shape[0]):
            var1 += measurement[i]*model[i]/(uncertainty[i]**2)
            var2 += (model[i]/uncertainty[i])**2
    return var1/var2

def giveLuminosity(I,I_err,lam,redshift,T_cold,T_hot,Ratio,Beta):
    # This method can use the FD (I), FD error(I_err) and redshift to determine the luminosity
    import numpy as np
    frequencies = np.linspace(3e11,3.75e13,(375-2))
    Lsun = 3.846e26
    intensities = np.zeros([len(frequencies)])
    model = tomModel((3.e8/lam),1,redshift,T_cold,T_hot,Ratio,Beta)
    amplitude = (1.e-26)*(frequencies[1]-frequencies[0])*giveAmplitude(model,I,I_err)*(4*np.pi*(Dlpen(redshift,giveAnswerInMeters=True)**2))/Lsun
    for i in range(len(frequencies)):
        intensities[i] = amplitude*tomModel(frequencies[i],1,redshift,T_cold,T_hot,Ratio,Beta)
    return intensities.sum()


##-----------------------------------------
## Give a cosmological distance
##-----------------------------------------

def Dlpen(redshift, giveAnswerInMeters = False):
	Mpc = 3.0857e22
	Dl = cosmo.luminosity_distance(redshift).value
	if giveAnswerInMeters:
		return Dl*Mpc
	else:
		return Dl

##-----------------------------------------
## Fitting a Gaussian
##-----------------------------------------

### Fitting Gaussians
def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


##-----------------------------------------
## Generate spectral lines
##-----------------------------------------

def LFIRtoSL(Luminosity,z,variance,COlines='Kamenetzky',lines='Bonato',giveNames='No'):
    '''
    Input: Luminosity [log (Lo)], redshift (z), variance (0 - no OR 1 - yes), COlines (Kamenetzky OR Rosenboom), lines (Bonato OR Spinoglio), giveNames ('Table', 'Tex' OR 'No').
    This will produce a list of frequencies [GHz] and flux densities [Jy], with or without a random variance for the SLs:
    Row 0 to 12: CO (1-0) to CO(13-12)
    Row 13 to 17: SIII, SiII, OIII, NIII, OI,
    Row 18 to 24: OIII, NII, OI, CII, CI, CI, NII
    giveNames toggles returning the names of the spectral lines in a list of strings.
    '''
    import numpy as np
    from scipy.stats import norm
    Dl = Dlpen(z,giveAnswerInMeters=True)
    Dlmpc = Dlpen(z,giveAnswerInMeters=False)
    c = 3.0e8
    velocity = 600. #km/s
    df = velocity/300000. #km/s
    lSun = 3.826e26
    outputArray = np.zeros([25,4])
    # Load spectral line intensities
    if COlines == 'Kamenetzky':
        slco = np.genfromtxt('./K17_Table7', skip_header=1, dtype=np.float, delimiter=", ", unpack = False)
    elif COlines == 'Rosenboom':
        slco = np.loadtxt('./COcoeff', dtype=np.float, delimiter=" ", unpack = False)
    else:
        print('Did not recognise the CO-lines library, will be using Kamenetzky')
        slco = np.genfromtxt('./K17_Table7', skip_header=1, dtype=np.float, delimiter=", ", unpack = False)
    if lines == 'Bonato':
        sl = np.loadtxt('./coeffBonato', dtype=np.float, delimiter="    ", unpack = False)
    elif lines == 'Spinoglio':
        sl = np.loadtxt('./coeff_spinoglio', dtype=np.float, delimiter=", ", unpack = False)
    else:
        print('Did not recognise the line-library, will be using Bonatos line estimates')
        sl = np.loadtxt('./coeffBonato', dtype=np.float, delimiter="    ", unpack = False)
    for i in range(13):
        outputArray[i,0] = ((1.e-9)*(i+1)*(115*(10**9))/(1+z)) # The CO lines 1-0 to 13-12
    outputArray[13,0] = ((1.e-9)*c/((1+z)*33.48e-6)) # SIII
    outputArray[14,0] = ((1.e-9)*c/((1+z)*34.82e-6)) # SiII
    outputArray[15,0] = ((1.e-9)*c/((1+z)*51.81e-6)) # OIII
    outputArray[16,0] = ((1.e-9)*c/((1+z)*57.32e-6)) # NIII
    outputArray[17,0] = ((1.e-9)*c/((1+z)*63.18e-6)) # OI
    outputArray[18,0] = ((1.e-9)*c/((1+z)*88.36e-6)) # OIII
    outputArray[19,0] = ((1.e-9)*c/((1+z)*121.9e-6)) # NII
    outputArray[20,0] = ((1.e-9)*c/((1+z)*145.5e-6)) # OI
    outputArray[21,0] = ((1.e-9)*c/((1+z)*157.7e-6)) # CII
    outputArray[22,0] = ((1.e-9)*c/((1+z)*370.5e-6)) # CI
    outputArray[23,0] = ((1.e-9)*c/((1+z)*609.6e-6)) # CI
    outputArray[24,0] = ((1.e-9)*c/((1+z)*205e-6))    # NII
    if variance != 0:
        # Create three random numbers per galaxy, used to create the log-normal distribution
        randvar = (np.random.random(3))
        # This creates a log-normal distribution, with a maximum deviation of 2 sigma
        randvar = norm.ppf(randvar*0.96 + 0.02)/norm.ppf(0.02)
        for i in range(13):
            outputArray[i,2] = (randvar[0])
        outputArray[13,2] = (randvar[2])
        outputArray[14,2] = (randvar[2])
        outputArray[15,2] = (randvar[0])
        outputArray[16,2] = (randvar[0])
        outputArray[17,2] = (randvar[1])
        outputArray[18,2] = (randvar[0])
        outputArray[19,2] = (randvar[0])
        outputArray[20,2] = (randvar[1])
        outputArray[21,2] = (randvar[1])
        outputArray[22,2] = (randvar[1])
        outputArray[23,2] = (randvar[1])
        outputArray[24,2] = (randvar[0])
        # Create three random numbers per galaxy, used to create the log-normal distribution
        randvar = (np.random.random(3))
        # This creates a log-normal distribution, with a maximum deviation of 2 sigma
        randvar = norm.ppf(randvar*0.96 + 0.02)/norm.ppf(0.02)
        for i in range(13):
            outputArray[i,3] = (randvar[0])
        outputArray[13,3] = (randvar[2])
        outputArray[14,3] = (randvar[2])
        outputArray[15,3] = (randvar[0])
        outputArray[16,3] = (randvar[0])
        outputArray[17,3] = (randvar[1])
        outputArray[18,3] = (randvar[0])
        outputArray[19,3] = (randvar[0])
        outputArray[20,3] = (randvar[1])
        outputArray[21,3] = (randvar[1])
        outputArray[22,3] = (randvar[1])
        outputArray[23,3] = (randvar[1])
        outputArray[24,3] = (randvar[0])
    if lines == 'Spinoglio':
        # This part calculates the lines according to the Spinoglio
        # Two parts of the variance are taken into account, one on the steepness of the L_FIR - L_SL relation (outputArray[i,3])
        # and one on the total amplitude of the L_FIR - L_SL relation
        for i in range(13,22):
            #outputArray = (sigma ^ variance) * (Lsun * (10^ L_Spin_var) / freq-width) * (10^26 [Jy] / 4 pi Dl^2)
            outputArray[i,1]= ((10**sl[i-13,3])**outputArray[i,2]) * (lSun*(10**((sl[i-13,0] + outputArray[i,3]*sl[i-13,1])* Luminosity - sl[i-13,2]))/(df*outputArray[i,0]*1.e9)) * ((10**26)/(Dl*Dl*4.*np.pi))
    else:
        # This part calculates the lines according to Bonato
        # This method only has the variance on the total amplitude, as the steepness is fixed to 1
        for i in range(13,25):
            #outputArray = variance * (lum / f-width) * (10^26 / 4 pi Dl^2)
            outputArray[i,1]= (sl[i-13,2]**outputArray[i,2]) * (lSun*(10**(sl[i-13,0] * Luminosity - sl[i-13,1]))/(df*outputArray[i,0]*1.e9)) * ((10**26)/(Dl*Dl*4.*np.pi))
    if COlines != 'Rosenboom':
        # This method relies on many galaxies, documented in Kamenetzky
        # They use the unwieldy units of Laccent, which with a bit of effort can give the S_CO
        # The amplitude and the steepness are fitted in the paper, and simulated here.
        for i in range(13):
            # outputArray = StaceyValue * mult.factor * df_CII / df_CO * flux_density_CII
            Laccent = (Luminosity - slco[i,2] - outputArray[i,3]*slco[i,3])/(slco[i,0] + outputArray[i,2]*slco[i,1])
            outputArray[i,1] = (10**Laccent) *((1+z)*((115.*(i+1))**2))/((3.25e7)*(velocity)*(Dlmpc*Dlmpc))
    else:
        # From Rosenboom's paper, we extract the ratios in luminosity of each line, and then relate that to the CII luminosity
        # The variation is exactly that of the CII-line. The CII / CO (1-0) relation is from Stacey.
        # There is no scatter inside the CO-ladder
        for i in range(13):
            outputArray[i,1] = (1./4100.)*slco[i] * (outputArray[21,0] / outputArray[i,0]) * outputArray[21,1]
    if lines == 'Spinoglio':
        i = 24
        # NII_205 -> I need more data on how to calculate this property
        outputArray[24,1] = (2**outputArray[i,2]) * (lSun*(10**(1.05 * Luminosity - 4.747))/(df*outputArray[i,0]*1.e9)) * ((10**26)/(Dl*Dl*4.*np.pi))
        #
        outputArray[22,1] = 0.4 * outputArray[3,1]
        outputArray[23,1] = 0.3 * outputArray[6,1]
    if giveNames == 'Table':
        listOfSLs = ['CO (1-0)','CO (2-1)','CO (3-2)','CO (4-3)','CO (5-4)','CO (6-5)','CO (7-6)','CO (8-7)','CO (9 - 8)','CO (10 - 9)','CO (11-10)','CO (12 - 11)','CO (13-12)', 'SIII 33', 'SiII 35', 'OIII 52', 'NIII 57', 'OI 63', 'OIII 88', 'NII 122', 'OI 145', 'CII 158', 'CI 370', 'CI 610', 'NII 205']
        return outputArray,listOfSLs
    elif giveNames =='Tex':
        listOfSLs = 'CO (1 - 0) & CO (2 - 1) & CO (3 - 2) & CO (4 - 3) & CO (5 - 4) & CO (6 - 5) & CO (7 - 6) & CO (8 - 7) & CO (9 - 8) & CO (10 - 9) & CO (11 - 10) & CO (12 - 11) & CO (13 - 12) & SIII 33 & SiII 35 & OIII 52 & NIII 57 & OI 63 & OIII 88 & NII 122 & OI 145 & CII 158 & CI 370 & CI 610 & NII 205'
        return outputArray,listOfSLs
    else:
        return outputArray

def spectrum(luminosity,redshift,fLow=220,fHigh=440,numFreqBins = 1500,linewidth=300):
    """
    Creates a galaxy spectrum consisting of a blackbody continuum and spectral
    lines. Keywords:
    luminosity: in units of log(L_fir [L_sol])
    redshift: z
    fLow, fHigh: minimum and maximum frequency that will be in the spectrum in units of GHz
    numFreqBins: amount of linearly spaced frequency bins at which the spectrum should be evaluated
    linewidth: width of spectral lines in units of km/s
    # Output: freqArray -> array frequencies in units of GHz
    # Output: spectrum -> array of the flux densities in the spectrum in units of Jy
    """
    margin = 10
    # Generate frequency array
    freqArray = np.linspace(fLow-margin, fHigh+margin, numFreqBins)
	# Create spectrum according to Bakx+2018
    spectrum = tomModel(freqArray*(1.e9),1,redshift,T_cold,T_hot,Ratio,Beta)
	# Normalize the flux to the given far-IR luminosity
    normLum = giveLuminosity(np.array([spectrum[0],spectrum[0]]),np.array([1,1]),((3.e8)/((1.e9)*np.array([freqArray[0],freqArray[0]]))),redshift,T_cold,T_hot,Ratio,Beta)
    spectrum = (10.**luminosity)*spectrum/normLum
	# Add the spectrum lines
    B,names = LFIRtoSL(luminosity,redshift,0,giveNames='Table')
    for i in range(len(B)):
        specLine = gaus(freqArray,B[i,1]*(600/linewidth),B[i,0],B[i,0]*linewidth/(3.e5))
        spectrum += specLine
    return freqArray,spectrum

def plotspectrum(freqArray, spectrumArray):
    """
    Plots the spectra generated by spectrum(). Keywords:
    freqArray: array of frequencies in GHz
    spectrumArray: array of flux densities at the corresponding frequencies in Jy
    """
    plt.figure()
    plt.plot(freqArray, spectrumArray)
    plt.grid()
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('Flux density [Jy]')
    plt.xlim(np.min(freqArray), np.max(freqArray))
    plt.ylim(bottom=0)
    plt.show()