"""
IRAC imaging of the cluster field taken before and after the SN discovery 
epoch

Nothing obviously there, though difference imaging doesn't work very well for
the rotated IRAC PSF

"""

import os
import glob

import numpy as np

import astropy.io.fits as pyfits 
import astropy.wcs as pywcs

AORS = ['r58289152', 'r58785536']

def sync():
    for aor in AORS:
        os.system(f'aws s3 sync s3://grizli-v1/IRAC/AORS/{aor}/ ./{aor}/')
        
def mosaic():
    
    from grizli import utils
    from drizzlepac.astrodrizzle import ablot 
    
    head, outputwcs = utils.make_wcsheader(ra=24.51568758, dec=-21.92554133, size=240, pixscale=1.0, get_hdu=False, theta=0)
    
    files = {}
    drz = {}
    
    for ch in ['ch1', 'ch2']:
        files[ch] = {}
        drz[ch] = {}
        
        for aor in AORS:
            files[ch][aor] = glob.glob(f'{aor}/{ch}/*/*fits.gz')
            files[ch][aor].sort()
            
            sci_list = [pyfits.open(file)['CBCD'].data for file in files[ch][aor]]
            unc_list = [pyfits.open(file)['CBUNC'].data for file in files[ch][aor]]
            N = len(sci_list)
            wht_list = [1/unc**2 for unc in unc_list]
            bcd_headers = [pyfits.open(file)['CBCD'].header for file in files[ch][aor]]

            headers = [pyfits.open(file)['WCS'].header for file in files[ch][aor]]
            wcs_list = [pywcs.WCS(head, relax=True) for head in headers]

            for i in range(N):
                mask = (unc_list[i] == 0) | (~np.isfinite(sci_list[i]))
                mask |= (~np.isfinite(unc_list[i]))
                
                sci_list[i] -= headers[i]['PEDESTAL']
                sci_list[i][mask] = 0
                wht_list[i][mask] = 0
                wcs_list[i].pscale = utils.get_wcs_pscale(wcs_list[i])
                
            _drz = utils.drizzle_array_groups(sci_list, wht_list, wcs_list, outputwcs=outputwcs, kernel='square', pixfrac=0., verbose=True)
            
            for iter in range(2):
                cr_wht = [wht*1 for wht in wht_list]
                # CR rejection
                for i in range(N):
                    _blt = ablot.do_blot(_drz[0], _drz[4], wcs_list[i], 1., coeffs=True, interp='poly5')
                    cr = (sci_list[i] - _blt)*np.sqrt(wht_list[i]) > 12
                    cr &= (_blt != 0)
                    cr_wht[i] *= (cr == 0)
                    
                _drz = utils.drizzle_array_groups(sci_list, cr_wht, wcs_list, outputwcs=outputwcs, kernel='square', pixfrac=0., verbose=True)
            
            for k in ['DATE_OBS', 'MJD_OBS', 'INSTRUME', 'CHNLNUM']:
                _drz[3][k] = bcd_headers[0][k]
                    
            pyfits.writeto(f'{aor}-{ch}_drz_sci.fits', data=_drz[0], header=_drz[3], overwrite=True)
            pyfits.writeto(f'{aor}-{ch}_drz_wht.fits', data=_drz[1], header=_drz[3], overwrite=True)
            drz[ch][aor] = _drz
            