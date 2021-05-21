import numpy as np
import sncosmo
import matplotlib.collections as mcoll
from astropy.table import Table
from matplotlib import pyplot as pl, cm, ticker

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)



__CONTOUR_RANGE__ = [[-2,2],[-1,8]]

class LensedSupernova( object ):
    def __init__( self, photfile, lensmodelfile,
                  tableformat='ascii.commented_header'):
        self.photometryfile = photfile
        self.lensmodelfile = lensmodelfile
        self.photometrytable = Table.read(
            self.photometryfile, format=tableformat)
        self.lensmodeltable = Table.read(
            self.lensmodelfile, format=tableformat)
        # TODO ?: read in each image as a separate sncosmo object?
        # self.snphot = sncosmo.read_lc(self.photometryfile)

    def _get_snphot_datum(self, snid, band, colname):
        snid = int(snid)
        band = band.lower()
        colname = colname.lower()
        idx = np.where( (self.photometrytable['id']==snid) &
                        (self.photometrytable['band']==band) )[0][0]
        return (self.photometrytable[colname][idx])


    def magnitude(self, snid, band):
        """Returns the measured AB magnitude in the specified band"""
        flux = self._get_snphot_datum(snid, band, 'flux')
        zp = self._get_snphot_datum(snid, band, 'zp')
        mag = -2.5*np.log10(flux) + zp
        return mag

    def magnitude_err_obs(self, snid, band):
        """Returns the observational magnitude error, based only
        on the statistical photometric uncertainty.
        (use magnitude_err_sys for the systematic error, incorporating
        lens modeling uncertainty)
        """
        flux = self._get_snphot_datum(snid, band, 'flux')
        fluxerr = self._get_snphot_datum(snid, band, 'fluxerr')
        magerr_stat = 1.0857 * fluxerr / flux
        return (magerr_stat)

    def color_and_err(self, snid, band1, band2):
        """Returns the observed color in magnitudes
        and (statistical only) observational error
        """
        mag1 = self.magnitude(snid, band1)
        mag2 = self.magnitude(snid, band2)
        flux1 = self._get_snphot_datum(snid, band1, 'flux')
        fluxerr1 = self._get_snphot_datum(snid, band1, 'fluxerr')
        flux2 = self._get_snphot_datum(snid, band2, 'flux')
        fluxerr2 = self._get_snphot_datum(snid, band2, 'fluxerr')
        fratioerr = np.sqrt(
            (fluxerr1/flux1)**2 + (fluxerr2/flux2)**2 )
        colorerr = 1.0857 * fratioerr
        return (mag1 - mag2, colorerr)

    def _get_lensmodel_datum(self, snid, lensmodelid, colname):
        snid = int(snid)
        lensmodelid = lensmodelid.upper()
        colname = colname.lower()
        idx = np.where( (self.lensmodeltable['modelid']==lensmodelid) &
                        (self.lensmodeltable['snid']==snid) )[0][0]
        return( self.lensmodeltable[colname][idx] )

    def magnification(self, snid, lensmodelid):
        """Returns the model magnification for the given SN ID
        and Lens model ID"""
        return (self._get_lensmodel_datum(snid, lensmodelid, 'mu'))

    def magnification_err(self, snid, lensmodelid):
        """Returns the lens model magnification uncertainty for the given SN ID
        and Lens model ID"""
        return (self._get_lensmodel_datum(snid, lensmodelid, 'mu_err'))

    def time_delay(self, snid, lensmodelid):
        """Returns the lens model magnification uncertainty for the given SN ID
        and Lens model ID"""
        return (self._get_lensmodel_datum(snid, lensmodelid, 'tdelay'))

    def model_corrected_magnitude_and_err(self, snid, lensmodelid, band):
        """Returns the magnification-corrected SN magnitude for the given
        SN image (snid), lens model (lensmodelid), and observed bandpass (band)
        as well as the composite uncertainty (quadratic sum of the
        observational uncertainty and lens model uncertainty in magnitudes)
        """
        mag = self.magnitude(snid, band)
        mag_err = self.magnitude_err_obs(snid, band)
        mu = self.magnification(snid, lensmodelid)
        mu_err = self.magnification_err(snid, lensmodelid)
        mag_corrected = mag + 2.5 * np.log10(mu)
        mag_corrected_err = np.sqrt( mag_err ** 2 + (2.5*np.log10(mu_err))**2 )
        return (mag_corrected, mag_corrected_err)

    def magnitude_err_sys(self, snid, band, bestlensmodel='E' ):
        """Returns the 'systematic uncertainty' estimate for the SN magnitude,
        due to uncertainty in the LENSTOOL model construction.  This is derived
        from the set of all lens model variants A-E.  The user-provided
        bestlensmodelid is taken as the 'fiducial' lens model, which the
        systematic uncertainty estimate is measured against
        """
        lensmodelidlist = np.unique(self.lensmodeltable['modelid'])
        deltamagnitudelist = []
        m_best, dm_best = self.model_corrected_magnitude_and_err(
            snid, bestlensmodel, band)
        for lensmodelid in lensmodelidlist:
            m, dm = self.model_corrected_magnitude_and_err(
                snid, lensmodelid, band)
            deltamagnitudelist.append(m-m_best)
        dm_sys_pos = np.max(deltamagnitudelist)
        dm_sys_neg = np.min(deltamagnitudelist)
        return ([[dm_sys_neg], [dm_sys_pos]])


    #@property
    #def tobs_range( self ):
    #    return( [ self.mjdobs-self.mjdpk-self.mjdpkerr,self.mjdobs-self.mjdpk+self.mjdpkerr ] )

    #@property
    #def tobs( self ):
    #    return( self.mjdobs-self.mjdpk )


#from pytools import colorpalette as cp
_COLOR1 = 'maroon' #cp.maroon
_COLOR2 = 'teal' #cp.teal
_COLOR3 = 'darkorange' #cp.darkorange

# Dictionary of sncosmo CCSN model names and their corresponding SN sub-type
SubClassDict_SNANA = {    'ii':{    'snana-2007ms':'IIP',  # sdss017458 (Ic in SNANA)
                                    'snana-2004hx':'IIP',  # sdss000018 PSNID
                                    'snana-2005gi':'IIP',  # sdss003818 PSNID
                                    'snana-2006gq':'IIP',  # sdss013376
                                    'snana-2006kn':'IIP',  # sdss014450
                                    'snana-2006jl':'IIP',  # sdss014599 PSNID
                                    'snana-2006iw':'IIP',  # sdss015031
                                    'snana-2006kv':'IIP',  # sdss015320
                                    'snana-2006ns':'IIP',  # sdss015339
                                    'snana-2007iz':'IIP',  # sdss017564
                                    'snana-2007nr':'IIP',  # sdss017862
                                    'snana-2007kw':'IIP',  # sdss018109
                                    'snana-2007ky':'IIP',  # sdss018297
                                    'snana-2007lj':'IIP',  # sdss018408
                                    'snana-2007lb':'IIP',  # sdss018441
                                    'snana-2007ll':'IIP',  # sdss018457
                                    'snana-2007nw':'IIP',  # sdss018590
                                    'snana-2007ld':'IIP',  # sdss018596
                                    'snana-2007md':'IIP',  # sdss018700
                                    'snana-2007lz':'IIP',  # sdss018713
                                    'snana-2007lx':'IIP',  # sdss018734
                                    'snana-2007og':'IIP',  # sdss018793
                                    'snana-2007ny':'IIP',  # sdss018834
                                    'snana-2007nv':'IIP',  # sdss018892
                                    'snana-2007pg':'IIP',  # sdss020038
                                    'snana-2006ez':'IIn',  # sdss012842
                                    'snana-2006ix':'IIn',  # sdss013449
                                },
                          'ibc':{    'snana-2004fe':'Ic',
                                     'snana-2004gq':'Ic',
                                     'snana-sdss004012':'Ic',  # no IAU ID
                                     'snana-2006fo':'Ic',      # sdss013195 PSNID
                                     'snana-sdss014475':'Ic',  # no IAU ID
                                     'snana-2006lc':'Ic',      # sdss015475
                                     'snana-04d1la':'Ic',
                                     'snana-04d4jv':'Ic',
                                     'snana-2004gv':'Ib',
                                     'snana-2006ep':'Ib',
                                     'snana-2007y':'Ib',
                                     'snana-2004ib':'Ib',   # sdss000020
                                     'snana-2005hm':'Ib',   # sdss002744 PSNID
                                     'snana-2006jo':'Ib',   # sdss014492 PSNID
                                     'snana-2007nc':'Ib',   # sdss019323
                                 },
                          'ia': {'salt2-extended':'Ia'},
                      }


SubClassDict_PSNID = {
           'ii':{ 's11-2004hx':'II','s11-2005lc':'IIP','s11-2005gi':'IIP','s11-2006jl':'IIP' },
           'ibc':{ 's11-2005hl':'Ib','s11-2005hm':'Ib','s11-2006fo':'Ic', 's11-2006jo':'Ib'},
           'ia': {'salt2-extended':'Ia'},
}


# Probability that a CC SN belongs to any given CC sub-class (Ib,Ic,IIP,IIL,IIn)
# from Li et al 2011a
ccSubClassProbs = {
           # For simulations containing only Type II sub-types
           'ii':{'IIP':0.7,'IIn':0.3 },
           # For simulations containing only Type Ib/c sub-types
           'ibc':{'Ib':0.46,'Ic':0.54},
           # For simulations containing all CC sub-types together
           'cc':{'Ib':0.19/0.76*0.46,'Ic':0.19/0.76*0.54,
                 'IIP':0.57/0.76*0.7,'IIn':0.57/0.76*0.3 },
           }


def mkobservationsTable( filterset='hst', orbits=10., tobs=0.  ):
    from simparam import gethstzp, gethstbgnoise
    import numpy as np
    from astropy.table import Table

    # medium band at peak observation set :
    if filterset == 'kfo':
        bandlist = ['f105w','f160w'] # for the lensed hst
    if filterset == 'hst' :
        bandlist = ['f350lp','f105w','f125w','f140w','f160w']
    elif filterset == 'jwst_z2' :
        bandlist = ['f140m','f162m','f182m','f210m',
                    'f150w','f200w',]
    elif filterset == 'jwst_z3' :
        bandlist = ['f140m','f162m','f182m','f210m','f250m','f300m','f335m','f360m',
                    'f150w','f200w','f277w','f356w']
    elif filterset == 'jwst_z4' :
        bandlist = ['f250m','f300m','f335m','f360m',
                    'f277w','f356w']
    elif filterset == 'jwst_z6' :
        bandlist = ['f250m','f300m','f335m','f360m',
                    'f277w','f356w' ]
    elif filterset == 'jwst_z8' :
        bandlist = ['f300m','f335m','f360m','f410m','f430m','f460m','f480m',
                    'f356w','f444w']

    if not np.iterable( orbits ) : orbits = np.ones(len(bandlist)) * orbits
    exptimelist = orbits * 2500

    timelist = np.zeros( len( bandlist ) ) + tobs
    zplist = [ gethstzp(band) for band in bandlist ]
    zpsyslist = [ 'ab' for i in range(len(bandlist)) ]
    # gainlist = [ gethstgain( et ) for et in exptimelist ]
    gainlist = [ et for et in exptimelist ]
    bgnoiselist = [ gethstbgnoise( band, et ) for band,et in zip(bandlist,exptimelist) ]
    # bgnoiselist = np.zeros( len(bandlist) )
    observations = Table({'band': bandlist, 'time': timelist, 'zp': zplist,
                          'zpsys': zpsyslist, 'gain': gainlist, 'exptime':exptimelist,
                          'skynoise': bgnoiselist,
                          })
    return( observations )


def mcsample( p, Ndraws, x0=None, mcsigma=0.05,
              Nburnin=100,  debug=False, *args, **kwargs ) :
    """ Crude metropolis-hastings monte carlo sampling funcion.

    The first argument is a callable function that defines
    the posterior probability at position x:  p(x).

    Positional arguments and optional keyword arguments for the function p
    may be provided at the end.  The function p will be called as
     p(x, *args, **kwargs).

    We construct a Markov Chain with  Ndraws  steps using the
    Metropolis-Hastings algorithm with a gaussian proposal distribution
    of stddev sigma.
    """
    from numpy import random
    if debug: import pdb; pdb.set_trace()

    # if user doesn't provide a starting point,
    # then draw an initial random position between 0 and 1
    if not x0 : x0 = random.uniform()
    xsamples = []
    istep = 0
    p0 = p(x0, *args, **kwargs)
    while len(xsamples) < Ndraws :
        # draw a new position from a Gaussian proposal dist'n
        x1 = random.normal( x0, mcsigma )
        p1 = p( x1, *args, **kwargs )
        # compare new against old position
        if p1>=p0 :
            # new position has higher probability, so
            # accept it unconditionally
            if istep>Nburnin : xsamples.append( x1 )
            p0=p1
            x0=x1
        else :
            # new position has lower probability, so
            # pick new or old based on relative probs.
            y = random.uniform( )
            if y<p1/p0 :
                if istep>Nburnin : xsamples.append( x1 )
                p0=p1
                x0=x1
            else :
                if istep>Nburnin : xsamples.append( x0 )
        istep +=1
    return( xsamples )

def pAv( Av, sigma=0, tau=0, R0=0, noNegativeAv=True ):
    """  Dust models:   P(Av)
    :param Av:
    :param sigma:
    :param tau:
    :param R0:
    :param noNegativeAv:
    :return:
    """
    import numpy as np
    if not np.iterable( Av ) : Av = np.array( [Av] )

    # gaussian core
    core = lambda sigma,av : np.exp( -av**2 / (2*sigma**2) )
    # Exponential tail
    tail = lambda tau,av : np.exp( -av/tau )

    if tau!=0 and noNegativeAv:
        tailOut = np.where( Av>=0, tail(tau,Av), 0 )
    elif tau!=0 :
        tailOut = tail(tau,Av)
    else :
        tailOut = np.zeros( len( Av ) )

    if sigma!=0 and noNegativeAv:
        coreOut = np.where( Av>=0, core(sigma,Av), 0 )
    elif sigma!=0 :
        coreOut = core(sigma,Av)
    else :
        coreOut = np.zeros( len( Av ) )

    if len(Av) == 1 :
        coreOut = coreOut[0]
        tailOut = tailOut[0]
    if sigma==0 : return( tailOut )
    elif tau==0 : return( coreOut )
    else : return( R0 * coreOut + tailOut )


class SncosmoSim( object ):
    """ An sncosmo SN population simulation,
    including parameters, cosmology, and light curves.
    """
    def __init__(self, sntype, observations=None, z_range=[1.8,2.2],
                 t0_range=[0,0], c_distrib=[0, 0.1], x1_distrib=[0, 1.],
                 nsim=100, perfect=True,
                 Om=0.3, H0=70, filterset='hst', mwdust=False, mwEBV=.018):
        """ Run a monte carlo sim using sncosmo to simulate <nsim> SNe
        of the given <sntype> over the given <z_range>.

        Simulates Type Ia SNe with the SALT2 model, and CC SNe with
        the SNANA CC templates.

        Observations are done at time t=0, unless specified otherwise in
        a user-defined observations table.

        Set perfect=True for noiseless "observations" of the simulated SNe.
        :return:
        """
        from astropy import cosmology
        import sncosmo
        from numpy.random import normal, uniform, choice
        import numpy as np

        self.sntype = sntype
        self.z_range = z_range
        self.nsim = nsim
        self.perfect = perfect

        if observations is None :
            observations = mkobservationsTable( filterset=filterset )
        self.observations = observations

        # Make a list of all the unique sncosmo source models available,
        # and assign a relative probability that any given simulated SN of this
        # type (CC or Ia) belongs to that subclass
        if sntype.lower() in ['cc','ii','ibc'] :
            subClassDict  = SubClassDict_SNANA[sntype.lower()]
            subClassProbs = ccSubClassProbs[sntype.lower()]
            #self.SourcenameSet = np.array( subClassDict.keys() )
            self.SourcenameSet = list(subClassDict.keys()) # kfo
            self.SubclassSet = np.array([ subClassDict[source] for source in self.SourcenameSet ])
            self.SubclassCount = np.array([ len(np.where(self.SubclassSet==subclass)[0])
                                         for subclass in self.SubclassSet ], dtype=float)
            self.SourceprobSet = np.array([ subClassProbs[subclass]
                                          for subclass in self.SubclassSet ]) / self.SubclassCount
            self.SourceprobSet /= self.SourceprobSet.sum()
        elif sntype.lower()=='ia' :
            # No sub-class divisions for SNIa
            self.SourcenameSet = np.array(['salt2-extended'])
            self.SubclassSet = np.array( ['Ia'] )
            self.SourceprobSet = np.array( [1] )
            self.SubclassCount = np.array( [1] )

        # load the O'Donnell 1994 dust model
        self.dust = sncosmo.OD94Dust()

        # Define an sncosmo SN model for each available source
        # kfo, mwdust defaults ~ False, will just use the host extinction effect  
        if mwdust:
            modelset = np.array([ sncosmo.Model(source=source, effects=[self.dust,self.dust],
                                    effect_names=['host','mw'], effect_frames=['rest','obs'])
                                 for source in self.SourcenameSet ])
        else:
            modelset = np.array([ sncosmo.Model(source=source, effects=[self.dust],
                                    effect_names=['host'], effect_frames=['rest'])
                                 for source in self.SourcenameSet ])

        # Define a cosmology
        # self.Om = Om
        # self.H0 = H0
        self.cosmo = cosmology.FlatLambdaCDM(Om0=Om, H0=H0)

        # For each simulated SN, draw random Av from distributions
        # as defined in Rodney et al 2014a :
        #   For SN Ia :  P(Av) = exp(-Av/0.33)
        #   For CC SN :  P(Av) = 4 * gauss(0.6) + exp(-Rv/1.7)
        if sntype=='Ia':
            tau,sigma,R0 = 0.33, 0, 0
        else :
            tau,sigma,R0 = 1.7, 0.6, 4
        self.Av = mcsample( pAv, nsim, tau=tau, sigma=sigma, R0=R0 )

        # For each simulated SN, draw a random Rv from a normal
        # distribution centered on 3.1 with width 0.5
        self.Rv = normal( 3.1, 0.5, nsim )
        self.Rv = np.where( self.Rv>0, self.Rv, 0.01 )

        # Convert Av and Rv to E(B-V) :
        # Rv = Av/EBV ==> EBV = Av/Rv
        self.EBV = self.Av / self.Rv

        # kfo convert, simplest assume mw Rv 3.1 each time, user provides ebv, solve extinction
        if mwdust:
            self.mwRv = 3.1 
            self.mwEBV = mwEBV
            self.mwAv = self.mwRv*self.mwEBV

        # TODO : draw the redshifts with a metropolis-hastings sampler to match
        #  a distribution defined based on the expected SN rate

        # Disabled : uniform redshift spacing
        # zlist = np.linspace( z_range[0], z_range[1], nsim )

        # Draw a random redshift from a uniform distribution
        self.z = uniform( low=z_range[0], high=z_range[1], size=nsim )

        lightcurvelist = []
        peakabsmagRlist = []
        modelparamlist = []
        subclasslist = []
        modelindexlist = []
        sourcenamelist = []
        t0list = []
        if sntype=='Ia':
            x0list = []
            x1list = []
            clist = []
        else :
            amplitudelist = []
        for isim in range(self.nsim):
            # Randomly draw an sncosmo model from the available list, according to
            # the predefined probability list, setting the SN sub-class for this
            # simulated SN
            imodel = choice( np.arange(len(modelset)), replace=True, p=self.SourceprobSet )
            model =  modelset[imodel]
            subclass = self.SubclassSet[imodel]

            z = self.z[isim]
            EBV = self.EBV[isim]
            Rv = self.Rv[isim]
            # kfo I am using Rv ~ 3.1; user provides optional color excess ebv 
            if mwdust:
                mwEBV = self.mwEBV
                mwRv = self.mwRv

            # Set the peak absolute magnitude according to the observed
            # luminosity functions, as defined in Table 3 of Graur:2014a;
            # and set the host extinction according to the 'mid' dust model
            # of Rodney:2014a.
            if subclass == 'Ia' :
                MR = normal( -19.37, 0.47 )
            elif subclass == 'Ib' :
                MR = normal( -17.90, 0.90 )
            elif subclass == 'Ic' :
                MR = normal( -18.30, 0.60 )
            elif subclass == 'IIP' :
                MR = normal( -16.56, 0.80 )
            elif subclass == 'IIL' :
                MR = normal( -17.66, 0.42 )
            elif subclass == 'IIn' :
                MR = normal( -18.25, 1.00 )
            model.set(z=z)
            model.set_source_peakabsmag( MR, 'bessellr', 'vega', cosmo=self.cosmo)

            modelindexlist.append( imodel )
            subclasslist.append( subclass )
            peakabsmagRlist.append( MR )
            sourcenamelist.append( list(self.SourcenameSet)[imodel] )
            if subclass =='Ia' :
                x0 = model.get('x0')
                # TODO : use bifurcated gaussians for more realistic x1,c dist'ns
                # kfo added x1_range and c_range arguments (default to commented)
                # x1 = normal(0., 1.)
                # c = normal(0., 0.1) 
                x1 = normal(x1_distrib[0], x1_distrib[1])
                c = normal(c_distrib[0], c_distrib[1])
                t0 = uniform( t0_range[0], t0_range[1] )
                # kfo 
                if mwdust:
                    modelparams = {'z':z, 't0':t0, 'x0':x0, 'x1':x1, 'c':c, 'hostebv':EBV, 'hostr_v':Rv,'mwebv':mwEBV,'mwr_v':mwRv}
                else:
                    modelparams = {'z':z, 't0':t0, 'x0':x0, 'x1':x1, 'c':c, 'hostebv':EBV, 'hostr_v':Rv}
                t0list.append( t0 )
                x0list.append( x0 )
                x1list.append( x1 )
                clist.append( c )
                t0list.append( t0 )
            else :
                amplitude = model.get('amplitude')
                t0 = uniform( t0_range[0], t0_range[1] )
                # kfo
                if mwdust:
                    modelparams = {'z':z, 't0':t0, 'amplitude':amplitude, 'hostebv':EBV, 'hostr_v':Rv,'mwebv':mwEBV,'mwr_v':mwRv }
                else:
                    modelparams = {'z':z, 't0':t0, 'amplitude':amplitude, 'hostebv':EBV, 'hostr_v':Rv }
                amplitudelist.append( amplitude )
                t0list.append( t0 )
            modelparamlist.append( modelparams )

            # Generate one simulated SN:
            snlc = sncosmo.realize_lcs(self.observations, model, [ modelparams ],
                                       thresh=None)#, perfect=perfect )
            lightcurvelist.append( snlc[0] )

        self.lightcurves = lightcurvelist
        self.t0 = np.array( t0list )
        self.modelindex = np.array( modelindexlist )
        self.sourcename = np.array( sourcenamelist )
        self.subclass = np.array( subclasslist )
        self.modelparam = np.array( modelparamlist )
        self.peakabsmagR = np.array( peakabsmagRlist )

        if sntype=='Ia':
            self.x0 = np.array( x0list )
            self.x1 = np.array( x1list )
            self.c  = np.array( clist )
        else :
            self.amplitude = np.array( amplitudelist )

        return


def scumsum( a ):
    """
    Sorted Cumulative Sum function :
    Construct an array "sumabove" such that the cell at index i in sumabove
    is equal to the sum of all cells from the input array "a" that have a
    cell value higher than a[i]
    """
    # Collapse the array into 1 dimension
    sumabove = a.ravel()

    # Sort the raveled array by descending cell value
    iravelsorted = sumabove.argsort( axis=0 )[::-1]

    # Reassign each cell to be the cumulative sum of all
    # input array cells with a higher value :
    sumabove[iravelsorted] = sumabove[iravelsorted].cumsum()

    # Now unravel back into shape of original array and return
    return( sumabove.reshape( a.shape ) )


def sncosmo_sim(snroot='nebra', z_range=[1.4,2.3], t0_range=[-20,20],
                c_distrib=[0.0, 0.1], x1_distrib=[0.5, 1.],
                filterset='hst', nsim=1000, verbose=True,
                clobber=False, mwdust=False, mwEBV=.018):
    """  Run sncosmo simulations for a color-color figure for SN Nebra
    """
    import os
    import _pickle as cPickle # kfo
    #import cPickle

    simIapkl='%s_SncosmoSim_Ia.pkl'%snroot
    simIIpkl='%s_SncosmoSim_II.pkl'%snroot
    simIbcpkl='%s_SncosmoSim_Ibc.pkl'%snroot

    if os.path.isfile( simIapkl ) and not clobber>1 :
        if verbose: print("Loading Ia simulation from pickle : %s"%simIapkl)
        fin = open( simIapkl, 'rb' )
        simIa = cPickle.load( fin )
        fin.close()
    else :
        if verbose: print("Running a new Ia simulation, then saving to pickle : %s"%simIapkl)
        simIa = SncosmoSim( 'Ia', z_range=z_range, t0_range=t0_range,
                            c_distrib=c_distrib, x1_distrib=x1_distrib, nsim=nsim, filterset=filterset, mwdust=mwdust, mwEBV=mwEBV)
        fout = open( simIapkl, 'wb' )
        cPickle.dump( simIa, fout, protocol=-1 )
        fout.close()

    if os.path.isfile( simIIpkl ) and not clobber>1 :
        if verbose: print("Loading II simulation from pickle : %s"%simIIpkl)
        fin = open( simIIpkl, 'rb' )
        simII = cPickle.load(fin)
        fin.close()
    else :
        if verbose: print("Running a new II simulation, then saving to pickle : %s"%simIIpkl)
        simII = SncosmoSim( 'II' , z_range=z_range, t0_range=t0_range, nsim=nsim, filterset=filterset,mwdust=mwdust,mwEBV=mwEBV )
        fout = open( simIIpkl, 'wb' )
        print('kfo cant pickle dict keys obj',simII)
        cPickle.dump( simII, fout, protocol=-1 )
        fout.close()

    if os.path.isfile( simIbcpkl ) and not clobber>1 :
        if verbose: print("Loading Ibc simulation from pickle : %s"%simIbcpkl)
        fin = open( simIbcpkl, 'rb' )
        simIbc = cPickle.load(fin)
        fin.close()
    else :
        if verbose: print("Running a new Ibc simulation, then saving to pickle : %s"%simIbcpkl)
        simIbc = SncosmoSim( 'Ibc' , z_range=z_range, t0_range=t0_range, nsim=nsim, filterset=filterset,mwdust=mwdust,mwEBV=mwEBV)
        fout = open( simIbcpkl, 'wb' )
        cPickle.dump( simIbc, fout, protocol=-1 )
        fout.close()

    return( simIa, simII, simIbc )

def _plot_colormag_singlesim(snsim, band1, band2, band3,
                               plotstyle='points',
                               nbins=None, colorcolor=True,contour_range=None,
                             ax=None, **plotargs):
    """ plot a color-mag diagram (kfo) x ~ (band1 - band2) y ~ (band3)
    :param snsim:
    :return:
    """
    if ax is None:
        ax = pl.gca()

    igood = np.where([np.all(snlc['flux']>0) for snlc in snsim.lightcurves])[0]
    # ibad = np.where([np.any(snlc['flux']<=0) for snlc in snsim.lightcurves])[0]
    lclist = [snsim.lightcurves[i] for i in igood]

    mag = np.array([
        -2.5*np.log10( np.ma.masked_less_equal(snlc['flux'],0,copy=False) )\
        + snlc['zp'] for snlc in lclist ])

    flt = np.array( [snlc['band'] for snlc in lclist] )
    i1 = np.where((flt == band1))
    i2 = np.where((flt == band2))
    i3 = np.where((flt == band3))

    # ax = pl.gca()
    # kfo, color mag x ~ (band1 - band2) y ~ (band3)
    if plotstyle=='points':
        plotargfinal = {'marker':'o', 'alpha':0.3, 'color':'darkorange', 'ls':' '}
        plotargfinal.update( **plotargs )
        ax.plot( mag[i1]-mag[i2], mag[i3], **plotargfinal )
    elif plotstyle.startswith('contour') or plotstyle=='gradient':
        xarray = mag[i2]-mag[i3]
        nsim = len(xarray)
        if nbins is None : nbins = int( np.sqrt( nsim  ) )
        if plotstyle.startswith('contour'):
            plotargfinal = {'levels':[0.0,0.68,0.95],'colors':['r','g','b'],
                            'ls':'-','alpha':0.5, 'extend':'neither'}
        else :
            plotargfinal = {'levels':np.arange(0.68,0.99,0.01),
                            'cmap':cm.Greys,'ls':'-', 'alpha':0.5,
                            'extend':'neither'}
        plotargfinal.update( **plotargs )

        # Plot filled contours, showing  the full extent of the population,
        # and contour lines containing 68% of the population.
        # First, bin the points into a 2-d histogram:
        # (Note that we reverse the x-y order here to get the binned arrays
        #  plotted in the correct direction )
        label = snsim.sntype # kfo gave SN type labels for a legend and added optional contour range. 
        count,y,x = np.histogram2d( mag[i3],mag[i1]-mag[i2],
                                    bins=nbins, range=contour_range) 

        # Renormalize relative to the sum of all SNe in this class :
        count /= count.sum()

        # Now set up an array 'cabove' such that  the cell value in cabove[i,j]
        # is equal to the sum of all cells that have a value higher than c[i,j]
        cabove = scumsum( count )

        # solid lines give probability contours at specified levels
        # (defaults to 0.68 for "1-sigma contours")
        #ax.contour( x[:-1], y[:-1], cabove, **plotargfinal )
        ax.contour( x[:-1], y[:-1], cabove, **plotargfinal )
        ax.contourf( x[:-1], y[:-1], cabove, **plotargfinal )
        plotargfinal['levels']=[0,.68]
        plotargfinal['alpha']=min(0.7,plotargfinal['alpha']*1.5)
        ax.contourf( x[:-1], y[:-1], cabove, **plotargfinal )

    ax.set_xlabel( '%s - %s' % (band1.upper(), band2.upper()))
    ax.set_ylabel( '%s' % (band3.upper()))
    #pl.legend()

    #ax.xaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    #ax.xaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    #ax.yaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    #ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    return( ax )

def _plot_colorcolor_singlesim(snsim, band1, band2, band3,
                               plotstyle='points',
                               nbins=None, colorcolor=True,**plotargs):
    """ plot a color-color diagram (or a color mag kfo)
    :param snsim:
    :return:
    """
    import numpy as np
    from matplotlib import pyplot as pl, cm, ticker

    igood = np.where([np.all(snlc['flux']>0) for snlc in snsim.lightcurves])[0]
    # ibad = np.where([np.any(snlc['flux']<=0) for snlc in snsim.lightcurves])[0]
    lclist = [snsim.lightcurves[i] for i in igood]

    mag = np.array([
        -2.5*np.log10( np.ma.masked_less_equal(snlc['flux'],0,copy=False) )\
        + snlc['zp'] for snlc in lclist ])

    flt = np.array( [snlc['band'] for snlc in lclist] )
    i1 = np.where((flt == band1))
    i2 = np.where((flt == band2))
    i3 = np.where((flt == band3))

    ax = pl.gca()


    if plotstyle=='points':
        plotargfinal = {'marker':'o', 'alpha':0.3, 'color':'darkorange', 'ls':' '}
        plotargfinal.update( **plotargs )
        ax.plot( mag[i1]-mag[i2], mag[i2]-mag[i3], **plotargfinal )
    elif plotstyle.startswith('contour') or plotstyle=='gradient':
        xarray = mag[i2]-mag[i3]
        nsim = len(xarray)
        if nbins is None : nbins = int( np.sqrt( nsim  ) )
        if plotstyle.startswith('contour'):
            plotargfinal = {'levels':[0.0,0.68,0.95],'colors':['r','g','b'],
                            'ls':'-','alpha':0.5, 'extend':'neither'}
        else :
            plotargfinal = {'levels':np.arange(0.68,0.99,0.01),
                            'cmap':cm.Greys,'ls':'-', 'alpha':0.5,
                            'extend':'neither'}
        plotargfinal.update( **plotargs )

        # Plot filled contours, showing  the full extent of the population,
        # and contour lines containing 68% of the population.
        # First, bin the points into a 2-d histogram:
        # (Note that we reverse the x-y order here to get the binned arrays
        #  plotted in the correct direction )
        count,y,x = np.histogram2d( mag[i2]-mag[i3],mag[i1]-mag[i2],
                                    bins=nbins, range=__CONTOUR_RANGE__ )

        # Renormalize relative to the sum of all SNe in this class :
        count /= count.sum()

        # Now set up an array 'cabove' such that  the cell value in cabove[i,j]
        # is equal to the sum of all cells that have a value higher than c[i,j]
        cabove = scumsum( count )

        # solid lines give probability contours at specified levels
        # (defaults to 0.68 for "1-sigma contours")
        #ax.contour( x[:-1], y[:-1], cabove, **plotargfinal )
        ax.contourf( x[:-1], y[:-1], cabove, **plotargfinal )

    ax.set_xlabel( '%s - %s' % (band1.upper(), band2.upper()))
    ax.set_ylabel( '%s - %s' % (band2.upper(), band3.upper()))

    #ax.xaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    #ax.xaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    #ax.yaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    #ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    return( ax )
    


def plotcontours( sim1, sim2, sim3=None, band1='f350lp', band2='f125w',
                  band3='f160w', nbins=None,colorcolor=True,contour_range=None,
                  ax=None, **plotargs ):
    """ Make a circle diagram, i.e. a med-wide band pseudo-color-color plot,
    showing both Type Ia and CC simulations over the given redshift range.
    :param snsim:
    :return:
    """
    if ax is None:
        ax = pl.gca()

    plotargs1 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR1], 'alpha':0.3 }
    plotargs1.update( **plotargs )

    plotargs2 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR2], 'alpha':0.3 }
    plotargs2.update( **plotargs )

    plotargs3 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR3], 'alpha':0.3 }
    plotargs3.update( **plotargs )
    # kfo
    if colorcolor:
        ax = _plot_colorcolor_singlesim(sim1, band1, band2, band3, nbins=nbins, plotstyle='contourf', ax=ax, **plotargs1)
        ax = _plot_colorcolor_singlesim(sim2, band1, band2, band3, nbins=nbins, plotstyle='contourf', ax=ax, **plotargs2)
        if sim3 is not None :
            ax = _plot_colorcolor_singlesim(sim3, band1, band2, band3, nbins=nbins, plotstyle='contourf', ax=ax, **plotargs3)
        return( ax )
    else: # kfo colormag, include optional contour_range 
        ax = _plot_colormag_singlesim(sim1, band1, band2, band3, nbins=nbins, plotstyle='contourf',contour_range=contour_range, ax=ax, **plotargs1)
        ax = _plot_colormag_singlesim(sim2, band1, band2, band3, nbins=nbins, plotstyle='contourf', contour_range=contour_range, ax=ax, **plotargs2)
        if sim3 is not None :
            ax = _plot_colormag_singlesim(sim3, band1, band2, band3, nbins=nbins, plotstyle='contourf', contour_range=contour_range, ax=ax, **plotargs3)
        return( ax )

def plotgradient( sim1, sim2, sim3=None, band1='f350lp', band2='f125w',
                  band3='f160w', nbins=None,colorcolor=True, **plotargs ):
    """ Make a circle diagram, i.e. a med-wide band pseudo-color-color plot,
    showing both Type Ia and CC simulations over the given redshift range.
    :param snsim:
    :return:
    """
    plotargs1 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR1], 'alpha':0.3 }
    plotargs1.update( **plotargs )

    plotargs2 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR2], 'alpha':0.3 }
    plotargs2.update( **plotargs )

    plotargs3 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR3], 'alpha':0.3 }
    plotargs3.update( **plotargs )
    # kfo
    if colorcolor:
        ax = _plot_colorcolor_singlesim(sim1, band1, band2, band3, nbins=nbins, plotstyle='gradient', **plotargs1)
        ax = _plot_colorcolor_singlesim(sim2, band1, band2, band3, nbins=nbins, plotstyle='gradient', **plotargs2)
        if sim3 is not None :
            ax = _plot_colorcolor_singlesim(sim3, band1, band2, band3, nbins=nbins, plotstyle='gradient', **plotargs3)
        return( ax )
    else: # kfo color mag
        ax = _plot_colormag_singlesim(sim1, band1, band2, band3, nbins=nbins, plotstyle='gradient', **plotargs1)
        ax = _plot_colormag_singlesim(sim2, band1, band2, band3, nbins=nbins, plotstyle='gradient', **plotargs2)
        if sim3 is not None :
            ax = _plot_colormag_singlesim(sim3, band1, band2, band3, nbins=nbins, plotstyle='gradient', **plotargs3)
        return( ax )




def plotpoints(sim1, sim2, sim3=None, band1='f350lp', band2='f125w',
               band3='f160w', colorcolor=True,**plotargs):
    """ Make a circle diagram, i.e. a med-wide band pseudo-color-color plot,
    showing both Type Ia and CC simulations over the given redshift range.
    :param snsim:
    :return:
    """
    # kfo
    if colorcolor:
        ax = _plot_colorcolor_singlesim(sim2, band1, band2, band3, plotstyle='points', marker='o', ls=' ', color=_COLOR2, **plotargs)
        if sim3 is not None :
            ax = _plot_colorcolor_singlesim(sim3, band1, band2, band3, plotstyle='points', marker='o', ls=' ', color=_COLOR3, **plotargs)
        ax = _plot_colorcolor_singlesim(sim1, band1, band2, band3, plotstyle='points', marker='o', ls=' ', color=_COLOR1, **plotargs)

        return( ax )
    else: # kfo color mag
        ax = _plot_colormag_singlesim(sim2, band1, band2, band3, plotstyle='points', marker='o', ls=' ', color=_COLOR2, **plotargs)
        if sim3 is not None :
            ax = _plot_colormag_singlesim(sim3, band1, band2, band3, plotstyle='points', marker='o', ls=' ', color=_COLOR3, **plotargs)
        ax = _plot_colormag_singlesim(sim1, band1, band2, band3, plotstyle='points', marker='o', ls=' ', color=_COLOR1, **plotargs)

        return( ax )



def colorline(x, y, do_times=False, t=None,
              t_step=20, z=None,
              cmap=pl.get_cmap('copper'),
              norm=pl.Normalize(0.0, 1.0),
              linewidth=3, alpha=1.0, ax=None):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """
    if ax is None:
        ax = pl.gca()

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)
    ax.add_collection(lc)

    if do_times:
        #tmp = list(zip([x,y,t])
        for tmp_t in np.arange(min(t),max(t),t_step):
            i = np.where(np.abs(t-tmp_t)==np.min(np.abs(t-tmp_t)))[0][0]
            x_loc,y_loc,t_val = x[i],y[i],t[i] # rf t_val at the color,mag value
            # none xycoords should default to same axis
            ax.annotate('t ~ {}'.format(int(t_val)),(x_loc,y_loc),fontsize=10)

    return lc


def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


# TODO: this wants some cleaning up
def get_sn_model_color_mag_time_data(
        source='salt2-extended', z=1.95,
        return_obsframe_time=True, do_dust=True,
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3),
        MR=-19.37, # Vega, Wang+2006 and Graur+2014; # alt normalization: MB=-19.25 (Richardson+ 2014)
        tmin=-20, tmax=100, tsteps=1000,
        param_dict=None,printt=False,ebv=.014,c=None):
    """ Generate the data for a SN color-magnitude vs time curve.
    Returns arrays of color, magnitude in F105W,  magnitude in F160W,
    and observer-frame time (days relative to peak brightness).
    """
    time=np.linspace(tmin, tmax, tsteps)

    # color curve m105w-m160w for a SN template
    if do_dust:
        dust = sncosmo.CCM89Dust()
        model = sncosmo.Model(source=source,
                              effects=[dust, dust],
                              effect_names=['host', 'mw'],
                              effect_frames=['rest', 'obs'])

        # look up ra dec if want to set mw dust
        # dustmap = sfdmap.SFDMap('/Users/kyleoconnor/Documents/GitHub/sfdmap/sfddata-master/')
        # ebv = dustmap.ebv(ra, dec)
        model.set(z=z,mwebv=ebv) # optionally include mwebv=ebv & **param_dict
    else:
        model = sncosmo.Model(source=source)
        model.set(z=z)
    if c != None: # color param in salt2, default is None, to show unc region need to provide arg
        model.set(c=c)

    if MR != None:
        model.set_source_peakabsmag( MR, 'bessellr', 'vega', cosmo=cosmo)
    #if MB != None:
    #    model.set_source_peakabsmag( MB, 'bessellb', 'vega', cosmo=cosmo)

    m160 = model.bandmag('f160w', 'ab', time)
    m105 = model.bandmag('f105w', 'ab', time)
    c = m105 - m160

    if printt:
        print(model.parameters,model.param_names)
    if return_obsframe_time is False:
        time /= (1+z)

    return [c,m105,m160,time]

def get_magnitude_arrays_from_sim(
        snsim, band1='f160w',band2='f105w',band3='f160w'):
    """Extract the magnitudes in specified bands from the sim object.
    Returns three arrays (possibly redundant), for each of the three bands.
    """
    igood = np.where([np.all(snlc['flux']>0) for snlc in snsim.lightcurves])[0]
    lclist = [snsim.lightcurves[i] for i in igood]
    mag = np.array([
        -2.5*np.log10( np.ma.masked_less_equal(snlc['flux'],0,copy=False) ) \
        + snlc['zp'] for snlc in lclist ])
    flt = np.array( [snlc['band'] for snlc in lclist] )
    i1 = np.where((flt == band1))
    i2 = np.where((flt == band2))
    i3 = np.where((flt == band3))
    m1 = mag[i1]
    m2 = mag[i2]
    m3 = mag[i3]
    return(m1,m2,m3)