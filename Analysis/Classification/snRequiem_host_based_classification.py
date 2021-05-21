import os
import numpy as np

from matplotlib import pyplot as plt

from astropy.io import fits
from astropy.table import QTable, Table
from astropy import units as u
from astropy import constants

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.05)

import sncosmo
from scipy.interpolate import interp1d

cgsfluxdensity = u.erg / u.Angstrom / u.cm**2 / u.s

_HOSTDATADIR_ =  "../Host/Photometry"
_HOSTPHOTFILE_ = "host_observed_total_mag_ab.dat"
_HOSTPHOTFORMAT_ = 'ascii.commented_header'
_BANDLIST_ = ['F225W', 'F275W', 'F336W', 'F390W', 'F555W', 'F606W','F775W',
              'UVF814W', 'F850LP', 'F105W','F110W','F125W','F140W', 'F160W',
              'IRAC1','IRAC2'] # no spitzer in sncosmo, wavelengths hardcoded below

_HOSTSEDFILE_  = "host_image2_sed.fits"
_HOSTSEDTEMPLATEFILE_ = "host_image2_template.fits"
_HOSTSEDFORMAT_ = 'fits'
_ZSRC_ = 1.95


class ObservedPhotometry():
    def __init__(self, datafile=_HOSTPHOTFILE_, datadir=_HOSTDATADIR_,
                 dataformat=_HOSTPHOTFORMAT_):
        """Observed photometry: AB mags"""
        # Read in the SN host galaxy total magnitudes (AB mag)
        self.filename = os.path.join(datadir, datafile)
        self.data = Table.read(self.filename, format=dataformat)

    def bandmag(self, image : int, bandpassname : str) ->float:
        """Return the measured AB mag and err for the given host galaxy
        image in the given bandpass"""
        imagename = f'Image_{image:d}'
        irow = self.data['image'].tolist().index(imagename)
        fcolname = f'f_{bandpassname.upper()}'
        ecolname = f'e_{bandpassname.upper()}'
        magab = self.data[fcolname][irow]
        magab_err = self.data[ecolname][irow]
        return(magab, magab_err)


    def bandflux(self, image : int, bandpassname : str) ->float:
        """Return the measured flux (fnu, Jansky) and uncertainty
        for the given host galaxy image in the given bandpass"""
        mag, magerr = self.bandmag(image, bandpassname)
        fnu = 3631 * 10**(-0.4 * mag)  * u.Jansky
        fnuerr =  0.92103 * magerr * fnu
        return( fnu, fnuerr )

    def observed_host_phot(self, image : int, bandlist=_BANDLIST_):
        """ Read in the SN host galaxy total magnitudes (AB mag) and the
        best-fit template SED (flambda in cgs units).
        Convert template SED to fnu in Jansky, and observed total AB mag to
        Jansky, filter by filter.
        Rescale the template SED to match the observed total AB mag.
        Apply a magnification correction, cosmological distance correction, and
        1+z redshift correction to get the rest-frame SED in Jansky.
        Multiply by the B and K bandpasses and integrate to get absolute mags in
        B and K.
        """
        # read in the best-fit template SED (flambda in cgs units).
        # Convert template SED to fnu in Jansky

        fluxlist = []
        fluxerrlist = []
        wavelist = []
        for bandpassname in bandlist:
            f, ferr = self.bandflux(image, bandpassname)
            if f > 1 * u.Jansky:
                continue
            fluxlist.append(f.value)
            fluxerrlist.append(ferr.value)
            wavelist.append(bandwave(bandpassname))

        return(np.array(wavelist)*u.Angstrom,
               np.array(fluxlist)*u.Jansky,
               np.array(fluxerrlist)*u.Jansky)


class ObservedSED():

    def __init__(self, datafile=_HOSTSEDFILE_, datadir=_HOSTDATADIR_,
                 format=_HOSTSEDFORMAT_):
        self.filename=os.path.join(datadir, datafile)
        self.data = QTable.read(self.filename)

        self.wave = self.data['filter_pivot']
        self.flambda = self.data['fobs'] #* u.erg / u.Angstrom / u.cm**2 / u.s
        self.fnu = jansky_from_flambda(self.flambda, self.wave)



class TemplateSED():

    def __init__(self, datafile=_HOSTSEDTEMPLATEFILE_, datadir=_HOSTDATADIR_,
                 format=_HOSTSEDFORMAT_):
        self.filename=os.path.join(datadir, datafile)
        self.data = QTable.read(self.filename)

        self.wave = self.data['template_wave']
        self.flambda = self.data['templf'] #* u.erg / u.Angstrom / u.cm**2 / u.s
        self.fnu = jansky_from_flambda(self.flambda, self.wave)

    def get_absolute_magnitude(self, bandpassname, zobs=_ZSRC_, mu=1.):
        """Returns the absolute magnitude in the given bandpass,
        assuming the observed SED is at the given redshift zobs and
        magnified by the given magnification mu."""

        # convenience shorthands for the filters we want
        if bandpassname=='B':
            bandpassname='bessellb'
        if bandpassname=='K':
            bandpassname='cspk'

        DL_Mpc = cosmo.luminosity_distance(zobs)

        # Apply the cosmological distance correction to rescale
        # the spectral flux density to 10 pc
        fnu_10pc = self.fnu * ( DL_Mpc.to(u.pc) / (10*u.pc) )**2

        # Apply a magnification correction
        fnu_10pc *= mu

        # Aply 1+z redshift correction to get the rest-frame wavelength.
        wave_restframe = self.wave / (1+zobs)

        # Get the absolute magnitudes in the given bandpass
        wave_eff = bandwave(bandpassname)
        interpolator = interp1d(wave_restframe, fnu_10pc)
        fnu_10pc_at_wave_eff = interpolator(wave_eff) * u.Jansky
        absmag = magnitude_from_fnu(fnu_10pc_at_wave_eff)

        return(absmag)


def bandwave(bandpassname : str) -> float:
    """Return the transmission-weighted mean wavelength of the
    given bandpass in Angstrom (via sncosmo except for IRAC)"""
    if bandpassname == 'IRAC1':
        return(35510)
    elif bandpassname == 'IRAC2':
        return(44930)
    bandpass = sncosmo.get_bandpass(bandpassname.upper())
    return(bandpass.wave_eff)

@u.quantity_input
def magnitude_from_fnu(fnu: u.jansky):
    """convert spectral flux density in Jansky to AB mag"""
    fnu_jansky = fnu.to(u.jansky)
    return -2.5*np.log10( fnu_jansky / (3631*u.Jansky) )

@u.quantity_input
def jansky_from_flambda(flambda: u.erg / u.Angstrom / u.cm**2 / u.s ,
                        wave : u.Angstrom):
    """Convert spectral flux density and wavelength from flambda in cgs
    to fnu in Jansky"""
    return (flambda * wave**2 / constants.c).to(u.Jansky)

@u.quantity_input
def magnitude_from_flambda(flambda: u.erg / u.Angstrom / u.cm**2 / u.s ,
                           wave : u.Angstrom):
    """Convert from spectral flux density flambda in cgs to AB mag"""
    fnu_jansky = ((3.34e4 * (wave.to(u.Angstrom)/u.Angstrom)**2
                   * flambda.to(cgsfluxdensity).value) ) * u.Jansky
    return -2.5*np.log10( fnu_jansky / (3631*u.Jansky) )

