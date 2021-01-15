#! /usr/bin/env python

#from stardust2.stardust import classify
import sys
sys.path.append('../../../starDust2')
import stardust
from stardust import classify
from astropy.table import Table
from copy import deepcopy
import numpy as np

import sncosmo
import pickle
from matplotlib import pyplot as plt

photdatatable = """#
# Final photometry and best available lens model (previously model H, now model E)
# flux is psf-fitting photometry, with associated uncertainty fluxerr
# 
# The magnification estimates `mu` have not been applied to the photometry
# the time delay estimates 'dt' have not been applied to the dates
#
# time   band    flux  fluxerr zp zpsys snid img   mu   muerr    dt
57588.03 f105w 0.0877  0.0165 23.9 ab     1  1.1  3.907 0.534   0.0
57588.03 f105w 2.3471  0.0160 23.9 ab     2  1.2  7.381 3.044 101.4
57588.03 f105w 0.1748  0.0157 23.9 ab     3  1.3  5.020 1.217  19.3
57587.97 f160w 0.6113  0.0433 23.9 ab     1  1.1  3.907 0.534   0.0
57587.97 f160w 3.5715  0.0453 23.9 ab     2  1.2  7.381 3.044 101.4
57587.97 f160w 1.1265  0.0442 23.9 ab     3  1.3  5.020 1.217  19.3
"""
snphotdata = Table.read(photdatatable, format='ascii.commented_header',
                        header_start=-1)


def get_lens_model_corrected_photometry(photdat):
    """ make a modified version of the photometry table,
    correcting for the lens model mu and dt"""

    photdatnew = deepcopy(photdat)

    mu = photdat['mu'].data
    muerr = photdat['muerr'].data
    dt = photdat['dt'].data
    f = photdat['flux'].data
    ferr = photdat['fluxerr'].data

    photdatnew['flux'] = f / mu
    photdatnew['fluxerr'] = (f / mu) * np.sqrt(
        (ferr / f) ** 2 + (muerr / mu) ** 2)
    photdatnew['time'] = photdat['time'] - dt

    return (photdatnew)

priors_from_host_data = {'Ia': 0.624,
                         'II': 0.358 * (57. / (57 + 19)),
                         'Ibc': 0.358 * (19. / (57 + 19))}

def run_stardust(verbose=False):
    photdatnew = get_lens_model_corrected_photometry(snphotdata)
    classification_results = classify.classify(
        photdatnew, zhost=1.95, zhosterr=0.001, t0_range=[57475.,57505.],
        zminmax=[1.94, 1.96], npoints=300, maxiter=5000, templateset='SNANA',
        nsteps_pdf=101, priors=priors_from_host_data,
        inflate_uncertainties=False, verbose=verbose)

    print("P(Ia)=%.4f"%classification_results['pIa'])
    print("P(II)=%.4f"%classification_results['pII'])
    print("P(Ib/c)=%.4f"%classification_results['pIbc'])

    pickle.dump( classification_results,
                 open( "snReqieum_stardust_classify_results.pkl", "wb" ) )

    print("Success. Pickled.")
    return(classification_results)


def run_stardust_quick(verbose=False):
    priors_from_host_data = {'Ia': 0.624,
                             'II': 0.358 * (57. / (57 + 19)),
                             'Ibc': 0.358 * (19. / (57 + 19))}

    snphotdata = Table.read(photdatatable,
                            format='ascii.commented_header',
                            header_start=-1)
    photdatnew = get_lens_model_corrected_photometry(snphotdata)
    classification_results = classify.classify(
        photdatnew, zhost=1.95, zhosterr=0.001, t0_range=[57475.,57505.],
        zminmax=[1.94,1.96], npoints=15, maxiter=50, templateset='PSNID',
        nsteps_pdf=26, priors=priors_from_host_data,
        inflate_uncertainties=False, verbose=verbose)

    print("P(Ia)=%.4f"%classification_results['pIa'])
    print("P(II)=%.4f"%classification_results['pII'])
    print("P(Ib/c)=%.4f"%classification_results['pIbc'])

    pickle.dump( classification_results,
                 open( "snReqieum_stardust_classify_results.pkl", "wb" ) )

    print("Success. Pickled.")
    return(classification_results)


def main():
    results = run_stardust(verbose=False)
    classify.plot_maxlike_fit(results['salt2-extended'], templateset='psnid')
    plt.savefig("/Users/rodney/Desktop/snRequiem_stardust_classify_Ia.png", overwrite=True)

    classify.plot_fits(results, nshow=3, templateset='psnid')
    plt.savefig("/Users/rodney/Desktop/snRequiem_stardust_classify_alltypes.png", overwrite=True)



