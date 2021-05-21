import numpy as np
import sncosmo
import datetime
from astropy.io import ascii

#https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/spectralresponse/
ch1 = ascii.read("CH1_IRAC.csv")
ch2 = ascii.read("CH2_IRAC.csv")

def IRACbands(plot=False):
    import sncosmo
    from astropy.io import ascii
    import astropy.units as u
    #https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/spectralresponse/
    # microns and electrons/photon
    ch1 = ascii.read("CH1_IRAC.csv")
    wavelength = ch1['col1'] # micron
    transmission = ch1['col2'] # electrons/photon
    band1 = sncosmo.Bandpass(wavelength, transmission, name='IRAC_Channel1',wave_unit=u.micron)
    sncosmo.register(band1,force=True)
    #print(band.minwave(),band.maxwave(),band.wave_eff)
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.rcParams.update({'font.size': 15})
        fig,ax = plt.subplots(figsize=(10.5,8.5))
        ax.plot(wavelength,transmission)
        ax.set_xlabel("Wavelength [micron]")
        ax.set_ylabel("Transmission [electrons/photon]")
        plt.savefig('transmission_IRAC_Ch1.pdf',bbox_inches='tight')
        
    ch2 = ascii.read("CH2_IRAC.csv")
    wavelength = ch2['col1'] # micron
    transmission = ch2['col2'] # electrons/photon
    band2 = sncosmo.Bandpass(wavelength, transmission, name='IRAC_Channel2',wave_unit=u.micron)
    sncosmo.register(band2,force=True)
    #print(band.minwave(),band.maxwave(),band.wave_eff)
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.rcParams.update({'font.size': 15})
        fig,ax = plt.subplots(figsize=(10.5,8.5))
        ax.plot(wavelength,transmission)
        ax.set_xlabel("Wavelength [micron]")
        ax.set_ylabel("Transmission [electrons/photon]")
        plt.savefig('transmission_IRAC_Ch2.pdf',bbox_inches='tight')
        
    return band1,band2


if __name__ == "__main__":
	IRACbands()
	# Ia model
	model=sncosmo.Model(source="salt2-extended")

	# Luminosity Functions: Goldstein 2019 Table 1, MB Vega, https://arxiv.org/pdf/1809.10147.pdf
	band,sys = 'bessellb','vega'
	MIa,sigmaIa = -19.23,0.1
	# Luminosity Functions: Converted to AB Mag System, http://www.astronomy.ohio-state.edu/~martini/usefuldata.html, Blanton et al. 2007
	# B-Band m_AB - m_Vega ~ -0.09
	band,sys = 'bessellb','ab'
	dm = -0.09
	MIa += dm
	# Consider the bright end for our purpose... can we get any constraint
	magIa = MIa - sigmaIa # bright

	# July-18/19-2016 all 3 images were observed within 5 arcsec of mrg0138 galaxy multiple images
	t_hst_images = datetime.date(2016,7,18)
	# Channel 1 and Channel 2 IR images from IRAC Spitzer on October-13-2016
	t_irac_images = datetime.date(2016,10,13)

	dt = t_irac_images - t_hst_images

	# Table 1 in mrg0138_supernova, lensing delays and magnifications
	# considering image 2
	# photometry that seems to be Ia age before peak
	age_phot = -20 # +20 or -11
	# lens model suggests magnification
	mu = 7 # pm 3 
	# lensing magnification
	magIa_lensed = magIa - 2.5*np.log10(mu)

	# Setting Ia at source redshift using lensing and luminosity function
	# No dust and No stretch/color at this point, just want to see if mag is near detectable 
	model.set(z=1.95, t0=0)
	model.set_source_peakabsmag(magIa_lensed,band,sys)

	print("Image 2 detected hst {}. Magnification mu ~ {}".format(t_hst_images,mu))
	print("Have model of Ia SN at z = 1.95, peak magnitude AB set as MIa = {} in {}".format(magIa,band))
	print("After magnification the peak magnitude MIa = {:.2f}".format(magIa_lensed))

	time = dt.days + age_phot 
	ch1_mag = model.bandmag('IRAC_Channel1','AB',time)
	ch2_mag = model.bandmag('IRAC_Channel2','AB',time)

	print("Image 2 at Ia age_hst ~ {} days (phot age observer frame relative to peak)".format(age_phot))
	print("IRAC CH1 and CH2 images on {} are {} days later than the hst detections on {}".format(t_irac_images,dt.days,t_hst_images))
	print("That puts Image 2 at Ia age_irac ~ {} days = (irac - hst) + age_hst".format(time))
	print("This is magnitudes AB of ch1 {:.2f} and ch2 {:.2f}".format(ch1_mag,ch2_mag))

	"""
	Archival IRAC images
	# https://sha.ipac.caltech.edu/applications/Spitzer/SHA/#id=SearchByRequestID&RequestClass=ServerRequest&DoSearch=true&SearchByRequestID.field.requestID=58785536,58289152&SearchByRequestID.field.includeSameConstraints=_none_&MoreOptions.field.prodtype=aor,pbcd&shortDesc=AORKEY&isBookmarkAble=true&isDrillDownRoot=true&isSearchResult=true
	# March-15-2016 --> 12 minute observation 30s frame time
	# Oct-13-2016 --> 18 minute observation  100s frame time

	Spitzer Performance Estimation Tools
	https://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/pet/pet.html
	Imaging Sensitivity https://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/pet/senspet/index.html
	# Using 30s frame Cryo IRAC and Warm IRAC
	"""
	# Point-source sensitivity micro-Jansky 1-sigma
	# 30s frame time ~ 30s exposure time/pixel
	PS_1sigma_warm_irac_ch1 = 2.22 # micro Jansky
	PS_1sigma_warm_irac_ch2 = 2.53 # micro Jansky
	PS_1sigma_cryo_irac_ch1 = 1.4 # micro Jansky
	PS_1sigma_cryo_irac_ch2 = 2.4 # micro Jansky

	print("Point Source 1sigma AB limiting mags")
	print("Channel Cryo ~ {:.2f}".format(-2.5*np.log10(PS_1sigma_cryo_irac_ch1*10**-6/3631)))
	print("Channel2 Cryo ~ {:.2f}".format(-2.5*np.log10(PS_1sigma_cryo_irac_ch2*10**-6/3631)))
	print("Channel Warm ~ {:.2f}".format(-2.5*np.log10(PS_1sigma_warm_irac_ch1*10**-6/3631)))
	print("Channel2 Warm ~ {:.2f}".format(-2.5*np.log10(PS_1sigma_warm_irac_ch2*10**-6/3631)))
