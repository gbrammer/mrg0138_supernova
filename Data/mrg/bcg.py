def bcg_surface_brightness():
    
    """
    Compute BCG surface brightness radial profiles
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    import astropy.units as u
    from astropy import constants as const
    from astropy.cosmology import WMAP9
    import astropy.io.fits as pyfits
    
    from astropy.modeling.models import Sersic1D, Gaussian1D, Polynomial1D
    from astropy.modeling.fitting import LevMarLSQFitter, SLSQPLSQFitter, SimplexLSQFitter
    
    im = pyfits.open('j013804m2156-f160w_drz_sci.fits')
    seg = pyfits.open('j013804m2156-ir_seg.fits.gz') 
    wcs = pywcs.WCS(im[0].header)
    
    z_bcg = 0.338
    dL = WMAP9.luminosity_distance(z_bcg).to(u.cm)
    
    # Pixel values to monochromatic luminosity per kpc2
    pix_kpc2 = (WMAP9.kpc_proper_per_arcmin(z_bcg).to(u.kpc/u.arcsec)*0.1*u.arcsec)**2
    to_Lsun = (im[0].header['PHOTFNU']*1e6*u.uJy)*4*np.pi*dL**2*(1+z_bcg)
    rest_wave = im[0].header['PHOTPLAM']/(1+z_bcg)
    nu = const.c/(rest_wave*u.Angstrom)
    to_Lsun *= nu
    # Units already converted
    to_Lsun = (to_Lsun.to(const.L_sun.unit)/const.L_sun).value
    to_Lsun_per_kpc2 = to_Lsun/pix_kpc2.value
    
    sh = im[0].data.shape
    yp, xp = np.indices(sh)
    
    # BCG coordinates
    ra, dec = 24.51569804, -21.92548062
    xc, yc = wcs.all_world2pix(np.array([[ra], [dec]]).T, 0)[0]
    seg_id = seg[0].data[int(yc), int(xc)]
    
    R = np.sqrt((xc-xp)**2+(yc-yp)**2)
    R_sn = np.sqrt(((xy_sn - np.array([xc, yc]))**2).sum(axis=1))

    theta = np.arctan2(xc-xp, yc-yp)
    theta_sn = np.arctan2(xc-xy_sn[:,0], yc-xy_sn[:,1])
    
    sn_coords = np.array([[24.51007534, -21.92734181],
                          [24.51231984, -21.92978758],
                          [24.51512533, -21.93066599]])
                          
    xy_sn = wcs.all_world2pix(sn_coords, 0)
    
    # By image
    img, pa, icl_sigma = 1, -1.23132, 200
    img, pa, icl_sigma = 2, -0.62954, 200
    img, pa, icl_sigma = 3, -0.09759, 200
    
    pa = theta_sn[img-1]
    colors = ['red', 'purple', 'orange']
    
    ang = 2/180*np.pi
    #dtheta = np.sin(ang)*R
    wedge = (np.abs(theta-pa) < ang) & (R < 300)
    wedge_seg = ((seg[0].data == 0) | (seg[0].data == seg_id))
    
    fitter = LevMarLSQFitter()
    r_eff, amp, n = 50, 0.612, 2
    prof = Sersic1D(r_eff=r_eff, amplitude=amp, n=n, fixed={'n':True}, bounds={'n':(0.5, 8)})
    
    if icl_sigma > 0:
        prof += Gaussian1D(mean=0, stddev=icl_sigma, fixed={'mean':True, 'stddev':False}, bounds={'amplitude':(0.001, 10), 'stddev':(50,800)})
        #prof += Polynomial1D(3, window=(0, R[wedge].max()))
        
    wedge_fit = wedge & wedge_seg #& (R < 50)
    _res = fitter(prof, R[wedge_fit], im[0].data[wedge_fit])
    
    _p = _res.parameters
    _n = _res.param_names
    _sersic = Sersic1D(amplitude=_p[_n.index('amplitude_0')], 
                       r_eff=_p[_n.index('r_eff_0')],
                       n=_p[_n.index('n_0')])

    _gau = Gaussian1D(amplitude=_p[_n.index('amplitude_1')], 
                                          mean=_p[_n.index('mean_1')],
                                          stddev=_p[_n.index('stddev_1')])

    # Make figure
    Rline = np.logspace(-2,3,1000)
        
    fig = plt.figure(figsize=[8, 4])
    ax = fig.add_subplot(121)
    
    ax.scatter(R[wedge & ~wedge_seg], im[0].data[wedge & ~wedge_seg]+dy, alpha=0.05, color='k', marker=',', s=5)
    ax.scatter(R[wedge & wedge_seg], im[0].data[wedge & wedge_seg]+dy, alpha=0.2, color='k', marker=',', s=5)

    pl = ax.plot(Rline, _res(Rline)+dy, linewidth=2, alpha=0.8, color=colors[img-1])
    ax.plot(Rline, _gau(Rline)+dy, linewidth=1, alpha=0.3, color=colors[img-1], zorder=100)
    ax.plot(Rline, _sersic(Rline)+dy, linewidth=1, alpha=0.3, color=colors[img-1], zorder=100)
    
    ax.loglog()
    ax.set_ylim(_res(10)*1.e-2, _res(10)*2+dy)  
    ax.set_ylabel('{1} image value + {0:.2f}, e-/s'.format(dy, im[0].header['FILTER']))
    #ax.set_yticklabels([0.1, 1])
    
    label = 'SN 1.{0}: \n'.format(img)
    label += '\n'+'$f$ = '+'{0:.3f} e/s/pix'.format(_res(R_sn[img-1]))
    label += '\n'+'log $\mu$ [$L_\odot/\mathrm{kpc}^2$]= '+'{0:.2f}'.format(np.log10(_res(R_sn[img-1])*to_Lsun_per_kpc2))
    
    ax.set_xlabel(r'$R_\mathrm{BCG}$, pix')
    
    ax.vlines(R_sn[img-1], ax.get_ylim()[0], ax.get_ylim()[1], color=pl[0].get_color(), linestyle='--', label=label)
    
    ax.set_xlim(10, R[wedge].max())
    ax.grid()
    ax.legend(loc='lower left')
    
    ww = wedge*1.
    ww[ww == 0] = np.nan

    ws = (wedge & wedge_seg)*1.
    ws[ws == 0] = np.nan
    
    ax = fig.add_subplot(122)
    
    slx, sly = slice(2282, 2695), slice(1467, 1889)  
    ax.imshow(im[0].data[sly, slx], vmin=-0.02, vmax=0.5)
    ax.imshow(ww[sly, slx], vmin=0, vmax=1, cmap='gray', alpha=0.1)
    ax.imshow(ws[sly, slx], vmin=0, vmax=1.5, cmap=colors[img-1].title()+'s', alpha=0.8)
    
    ax.scatter(xy_sn[img-1,0]-slx.start, xy_sn[img-1,1]-sly.start, marker='o', fc='None', ec=colors[img-1], s=100)

    ax.scatter(xc-slx.start, yc-sly.start, marker='x', color='k', s=20)
    
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    fig.tight_layout(pad=0.5)
    
    fig.savefig('sn_image_{0}_mu.png'.format(img))
