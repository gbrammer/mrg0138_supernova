"""
Discovery figure, Fig. 1
"""
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def for_present():

    from grizli.pipeline import auto_script
    from mastquery import overlaps
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
    import numpy as np      
    
    root = 'j013804m2156'

    x0, y0= 2465, 1740
    
    NY = 400
    NX = 400

    slx = slice(x0-NX, x0+NX)
    sly = slice(y0-NY, y0+NY)
    
    mi = -0.01
    
    si = 4.8*NX/200
    
    _rgb = auto_script.field_rgb(root=root, filters=['f105w','f125w','f160w'], HOME_PATH=None, xyslice=(slx, sly), show_ir=False, add_labels=False, output_format='png', suffix='.xon', rgb_min=mi, tick_interval=10, xsize=si, output_dpi=int(2*NX/si), tickparams={'axis': 'both', 'colors': 'w', 'which': 'both', 'length':0})
    
    rgb_scl = np.array([1.06,1.03,1.03])/1.03
    
    _rgb = auto_script.field_rgb(root=root, filters=['f110w','f125w','f140w'], HOME_PATH=None, xyslice=(slx, sly), show_ir=False, add_labels=False, output_format='png', rgb_scl=rgb_scl, suffix='.xoff', rgb_min=mi, tick_interval=10, xsize=si, output_dpi=int(2*NX/si), tickparams={'axis': 'both', 'colors': 'w', 'which': 'both', 'length':0})
    
    
def rgb():
    
    from grizli.pipeline import auto_script
    from mastquery import overlaps
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
    import numpy as np      
    
    from grizli import utils
<<<<<<< HEAD
=======

    plt.rcParams.update({'font.size':10})
>>>>>>> 8dafb3ee4ba770e7c3cd5c2541a900edcc2f58e2
    
    root = 'j013804m2156'

    x0, y0= 2465, 1740
    
    NY = 200
    NX = 200

    slx = slice(x0-NX, x0+NX)
    sly = slice(y0-NY, y0+NY)
    
    mi = -0.01
    
    si = 4.8
    
    _rgb = auto_script.field_rgb(root=root, filters=['f105w','f125w','f160w'], HOME_PATH=None, xyslice=(slx, sly), show_ir=False, add_labels=False, output_format='png', suffix='.on', rgb_min=mi, tick_interval=10, xsize=si, output_dpi=int(2*NX/si))

    _rgb = auto_script.field_rgb(root=root, filters=['f110w','f125w','f140w'], HOME_PATH=None, xyslice=(slx, sly), show_ir=False, add_labels=False, output_format='png', rgb_scl=[1.06,1,1.02], suffix='.off', rgb_min=mi, tick_interval=10, xsize=si, output_dpi=int(2*NX/si))
    
    im = pyfits.open(f'{root}-f160w_drz_sci.fits')
    wcs = pywcs.WCS(im[0].header, relax=True)
    wsl = wcs.slice((sly, slx))
    
    coo = utils.read_catalog('coords.txt')
    N, inter = 10, 0.5
    
    N, inter = 20, 1.0

    xy0 = np.array(wcs.all_world2pix(coo['ra'], coo['dec'], 0)).T
    
    xy = np.cast[int](np.round(wcs.all_world2pix(coo['ra'], coo['dec'], 0))).T
    
    for i in [4,5,6,7]:
        
            xi, yi = xy[i]
            slx = slice(xi-N, xi+N)
            sly = slice(yi-N, yi+N)
        
            _rgb = auto_script.field_rgb(root=root, filters=['f105w','f125w','f160w'], HOME_PATH=None, xyslice=(slx, sly), show_ir=False, add_labels=False, output_format='png', suffix='.on{0}'.format(coo['label'][i]), rgb_min=mi, tick_interval=inter, xsize=si/3.)

            _rgb = auto_script.field_rgb(root=root, filters=['f110w','f125w','f140w'], HOME_PATH=None, xyslice=(slx, sly), show_ir=False, add_labels=False, output_format='png', rgb_scl=[1.06,1,1.02], suffix='.off{0}'.format(coo['label'][i]), rgb_min=mi, tick_interval=inter, xsize=si/3.)
          
    
    # Label offsets
    #      H1  H2 H3   H4 SN3 SN2  SN1  SN4   SN5 3.1 3.2 3.3 3.4
    dx = [-30, 10, 0, -28,  0, 5,  -28, -27] + [0, 0, 0, -18, 0]
    dy = [-8, 14, 10,  0,  10, 10, -15,  10] + [15, 15, 15, 0, 0]
        
    ny = 3
    ny = 4
    
    if ny == 3:
        fig = plt.figure(figsize=(8,3/5.*8))
    else:
        fig = plt.figure(figsize=(8,4/6.*8))
        
    gs = GridSpec(ny, 2+ny, figure=fig)
    
    ax = fig.add_subplot(gs[:,:ny])
    #ax.plot([0,1])
    img = plt.imread('j013804m2156.on.png')
    sh = img.shape
    xyi = (xy - np.array([x0-NX, y0-NY]))/np.array([NX, NY])/2*np.array(sh[:-1])
    xyi0 = (xy0 - np.array([x0-NX, y0-NY]))/np.array([NX, NY])/2*np.array(sh[:-1])
    #dy = 20
    pixgrow = sh[0]/2/NX
    
    ax.imshow(img, origin='upper')
    ax.axis('off')
    ax.text(0.03, 0.97, 'a)', ha='left', va='top', transform=ax.transAxes, color='w')
    
    colors = ['w']*7 + ['pink'] + ['paleturquoise'] + ['w']*4

    for j, i in enumerate(range(8+5)):#[4,5,6]):            
        # if i == 7:
        #     col = 'pink'
        # elif i > 7:
        #     col = '0.5'
        # else:
        #     col = 'w'
        
        if i > 9:
            fs = 10
        else:
            fs = 11
            
        ax.text(xyi[i][0]+dx[i], sh[0]-xyi[i][1]+dy[i], coo['label'][i], ha='center', va='top',
                color=colors[i], fontsize=fs)
    
    labels = 'bcdefghijk'
    
    iters = [6,5,4]
    if ny == 4:
        iters += [7]
        
    for j, i in enumerate(iters):

        if i == 7 :
            col = 'pink'
        elif i == 8:
            col = 'paleturquoise'

        else:
            col = 'w'

        ax = fig.add_subplot(gs[j,ny])
        img = plt.imread('{0}.on{1}.png'.format(root, coo['label'][i]))
        ax.imshow(img, origin='upper')
        ax.axis('off')
        
        ax.text(0.05, 0.95, '{0}): {1}'.format(labels[j], coo['label'][i]), 
                ha='left', va='top', color=colors[i], transform=ax.transAxes)
        
        di = 0.02
        if j == 2:
            ax.text(0.06+di, 0.07, r'$y_\mathrm{{105}}$', ha='left', va='bottom', color='w', transform=ax.transAxes)
            ax.text(0.30+di, 0.07, r'$J_\mathrm{{125}}$', ha='left', va='bottom', color='0.7', transform=ax.transAxes)
            ax.text(0.49+di, 0.07, r'$H_\mathrm{{160}}$', ha='left', va='bottom', color='w', transform=ax.transAxes)
            
            ax.text(0.06+di, 0.18, '2016', ha='left', va='bottom', color='w',
                    transform=ax.transAxes)
            
        ax = fig.add_subplot(gs[j,ny+1])
        img = plt.imread('{0}.off{1}.png'.format(root, coo['label'][i]))
        ax.imshow(img, origin='upper')
        ax.axis('off')
        
        ax.text(0.05, 0.95, '{0})'.format(labels[j+ny], coo['label'][i]), 
                ha='left', va='top', color=colors[i], transform=ax.transAxes)
        
        if j == 2:
            ax.text(0.06+di, 0.07, r'$y_\mathrm{{110}}$', ha='left', va='bottom', color='w', transform=ax.transAxes)
            ax.text(0.3+di, 0.07, r'$J_\mathrm{{125}}$', ha='left', va='bottom', color='w', transform=ax.transAxes)
            ax.text(0.49+di, 0.07, r'$H_\mathrm{{140}}$', ha='left', va='bottom', color='w', transform=ax.transAxes)
            
            ax.text(0.06+di, 0.18, '2019', ha='left', va='bottom', color='w',
                    transform=ax.transAxes)
    
    # Error ellipse for SN4
    if 1:
        thet = np.linspace(0, 2*np.pi, 256)
        xt = np.cos(thet)*4.95*pixgrow
        yt = np.sin(thet)*1.14*pixgrow
        phi = (-126.5)/180*np.pi
        _mat = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
        xyreg = np.dot(np.array([xt, yt]).T, _mat)
        xreg = xyreg[:,0] + xyi0[7,0]
        yreg = sh[0]-(xyreg[:,1] + xyi0[7,1])
        fig.axes[0].plot(xreg, yreg, color=colors[7])
        
        # Subregion
        pixgrow2 = img.shape[0]/2/N
        xt = np.cos(thet)*4.95
        yt = np.sin(thet)*1.14
        phi = (-126.5)/180*np.pi
        _mat = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
        xyreg = np.dot(np.array([xt, yt]).T, _mat)
        xreg = (xyreg[:,0] + xy0[7,0] - xy[7][0] + N)*pixgrow2
        yreg = img.shape[0] - (xyreg[:,1] + xy0[7,1] - xy[7][1] + N)*pixgrow2
        
        for a in fig.axes[-2:]:
            a.plot(xreg, yreg, color=colors[7])

    # x for SN5
    if 1:
        fig.axes[0].plot(xyi0[8,0], sh[0] - xyi0[8,1], marker='x', ms=6, color=colors[8])    
        
    # Error ellipse for SN5
    if 0:
        thet = np.linspace(0, 2*np.pi, 256)
        xt = np.cos(thet)*0.35*pixgrow
        yt = np.sin(thet)*0.05*pixgrow
        phi = (-126.41)/180*np.pi
        _mat = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
        xyreg = np.dot(np.array([xt, yt]).T, _mat)
        xreg = xyreg[:,0] + xyi0[8,0]
        yreg = sh[0]-(xyreg[:,1] + xyi0[8,1])
        fig.axes[0].plot(xreg, yreg, color=colors[8])      

            
    gs.tight_layout(fig, pad=0.1)
    fig.savefig('xfig1_layout.pdf', dpi=150)   
    plt.close('all')
        
def phot_limits():
    """
    Compute photometric limits in circular apertures
    """
    import sep
    
    rd = np.array([(24.51007534, -21.92734181),
          (24.51231984, -21.92978758),
          (24.51512533, -21.93066599)])
    
    filts = {}
    for filt in ['f110w','f125w','f140w','f160w']:
        sci = pyfits.open(f'j013804m2156-{filt}_drz_sci.fits')
        
        wht = pyfits.open(f'j013804m2156-{filt}_drz_wht.fits')
        wcs = pywcs.WCS(sci[0].header)
        
        psf = pyfits.open(f'j013804m2156-{filt}_psf.fits')
        
        filts[filt] = {'im':sci, 'sci':sci[0].data.byteswap().newbyteorder(), 'wht':wht[0].data.byteswap().newbyteorder(), 'wcs':wcs, 'psf':psf['PSF','DRIZ1'].data.byteswap().newbyteorder()}
        
    # Bottom image
    
    xpi = filts['f140w']['wcs'].all_world2pix(rd, 0)
    xi = np.cast[int](np.round(xpi))
    N = 10
        
    filt = 'f140w'
    i = 2
    xi[2][1] -= 4
    xi[0][0] -= 4
    
    for filt in ['f110w','f125w','f140w']:
        slx = slice(xi[i][0]-N, xi[i][0]+N)
        sly = slice(xi[i][1]-N, xi[i][1]+N)
    
        subsci = filts[filt]['sci'][sly, slx]
    
        yprof = subsci.mean(axis=1)
        xprof = (subsci.T-yprof).T.mean(axis=0)
    
        yp, xp = np.indices(subsci.shape)-N
        _A = [subsci.flatten()*0.+1]
        order = 5
        for j in range(1,order+1):
            if i == 2:
                _A.append(yp.flatten()**j)
            else:
                _A.append(xp.flatten()**j)
    
        _A = np.array(_A).T
        _c = np.linalg.lstsq(_A, subsci.flatten(), rcond=-1)
        _m = _A.dot(_c[0]).reshape(subsci.shape)
    
        var = 1/filts[filt]['wht'][sly, slx]
        xc = xpi-xi+N

        ZP, ee = get_hst_ee(filt) 
        to_ujy = 10**(-0.4*(ZP-23.9))
        
        _phot = sep.sum_circle((subsci - _m).astype(np.float32)*to_ujy, [xc[i,0]], [xc[i,1]], [3.5], subpix=0, var=var*to_ujy**2)
        
        ZP = 23.9
        
        print('{6:.2f} {0} {1:4.2f}±{2:.2f}  {3:.2f} {4:.2f} {5:.2f}'.format(filt, _phot[0][0], _phot[1][0], ZP - 2.5*np.log10(_phot[0][0]/ee), 2.5/np.log(10)*_phot[1][0]/_phot[0][0], ZP - 2.5*np.log10(_phot[1][0]/ee*5), filts[filt]['im'][0].header['EXPSTART']))
    
        yp, xp = np.indices(subsci.shape)-N
        _A = [subsci.flatten()*0.+1]
        order = 5
        for j in range(1,order+1):
            if i == 2:
                _A.append(yp.flatten()**j)
            else:
                _A.append(xp.flatten()**j)
    
        _A = np.array(_A).T
        _c = np.linalg.lstsq(_A, subsci.flatten(), rcond=-1)
        _m = _A.dot(_c[0]).reshape(subsci.shape)
    
        var = 1/filts[filt]['wht'][sly, slx]
        xc = xpi-xi+N

        ZP, ee = get_hst_ee(filt) 
        to_ujy = 10**(-0.4*(ZP-23.9))
        
        _phot = sep.sum_circle((subsci - _m).astype(np.float32)*to_ujy, [xc[i,0]], [xc[i,1]], [3.5], subpix=0, var=var*to_ujy**2)
        
        ZP = 23.9
        
        print('{6:.2f} {0} {1:4.2f}±{2:.2f}  {3:.2f} {4:.2f} {5:.2f}'.format(filt, _phot[0][0], _phot[1][0], ZP - 2.5*np.log10(_phot[0][0]/ee), 2.5/np.log(10)*_phot[1][0]/_phot[0][0], ZP - 2.5*np.log10(_phot[1][0]/ee*5), filts[filt]['im'][0].header['EXPSTART']))
    
