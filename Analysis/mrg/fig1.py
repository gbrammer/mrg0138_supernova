"""
Discovery figure, Fig. 1
"""

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
    
    xy = np.cast[int](np.round(wcs.all_world2pix(coo['ra'], coo['dec'], 0))).T
    
    for i in [4,5,6]:
        
        xi, yi = xy[i]
        slx = slice(xi-N, xi+N)
        sly = slice(yi-N, yi+N)
        
        _rgb = auto_script.field_rgb(root=root, filters=['f105w','f125w','f160w'], HOME_PATH=None, xyslice=(slx, sly), show_ir=False, add_labels=False, output_format='png', suffix='.on{0}'.format(coo['label'][i]), rgb_min=mi, tick_interval=inter, xsize=si/3.)

        _rgb = auto_script.field_rgb(root=root, filters=['f110w','f125w','f140w'], HOME_PATH=None, xyslice=(slx, sly), show_ir=False, add_labels=False, output_format='png', rgb_scl=[1.06,1,1.02], suffix='.off{0}'.format(coo['label'][i]), rgb_min=mi, tick_interval=inter, xsize=si/3.)
        
    
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    
    # Label offsets
    dx = [-30, 10, 0, -10, 0, 0, -38]
    dy = [-8, 14, 10, 5, 10, 10, -15]
    
    fig = plt.figure(figsize=(8,3./5*8))
    gs = GridSpec(3, 5, figure=fig)
    
    ax = fig.add_subplot(gs[:,:3])
    #ax.plot([0,1])
    img = plt.imread('j013804m2156.on.png')
    sh = img.shape
    xyi = (xy - np.array([x0-NX, y0-NY]))/np.array([NX, NY])/2*np.array(sh[:-1])
    #dy = 20
    
    ax.imshow(img, origin='upper')
    ax.axis('off')
    ax.text(0.03, 0.97, 'a)', ha='left', va='top', transform=ax.transAxes, color='w')
    
    for j, i in enumerate(range(7)):#[4,5,6]):            
        ax.text(xyi[i][0]+dx[i], sh[0]-xyi[i][1]+dy[i], coo['label'][i], ha='center', va='top', color='w', fontsize=8)
    
    labels = 'bcdefg'
    
    for j, i in enumerate([6,5,4]):
        ax = fig.add_subplot(gs[j,3])
        img = plt.imread('{0}.on{1}.png'.format(root, coo['label'][i]))
        ax.imshow(img, origin='upper')
        ax.axis('off')
        
        ax.text(0.05, 0.95, '{0}): {1}'.format(labels[j], coo['label'][i]), 
                ha='left', va='top', color='w', transform=ax.transAxes)
        
        di = 0.02
        if j == 2:
            ax.text(0.06+di, 0.07, r'$y_\mathrm{{105}}$', ha='left', va='bottom', color='w', transform=ax.transAxes)
            ax.text(0.28+di, 0.07, r'$J_\mathrm{{125}}$', ha='left', va='bottom', color='0.7', transform=ax.transAxes)
            ax.text(0.49+di, 0.07, r'$H_\mathrm{{160}}$', ha='left', va='bottom', color='w', transform=ax.transAxes)
            
            ax.text(0.06+di, 0.18, '2016', ha='left', va='bottom', color='w',
                    transform=ax.transAxes)
            
        ax = fig.add_subplot(gs[j,4])
        img = plt.imread('{0}.off{1}.png'.format(root, coo['label'][i]))
        ax.imshow(img, origin='upper')
        ax.axis('off')
        
        ax.text(0.05, 0.95, '{0})'.format(labels[j+3], coo['label'][i]), 
                ha='left', va='top', color='w', transform=ax.transAxes)
        
        if j == 2:
            ax.text(0.06+di, 0.07, r'$y_\mathrm{{110}}$', ha='left', va='bottom', color='w', transform=ax.transAxes)
            ax.text(0.28+di, 0.07, r'$J_\mathrm{{125}}$', ha='left', va='bottom', color='w', transform=ax.transAxes)
            ax.text(0.49+di, 0.07, r'$H_\mathrm{{140}}$', ha='left', va='bottom', color='w', transform=ax.transAxes)
            
            ax.text(0.06+di, 0.18, '2019', ha='left', va='bottom', color='w',
                    transform=ax.transAxes)
                            
    gs.tight_layout(fig, pad=0.1)
    fig.savefig('fig1_layout.pdf', dpi=150)   
    plt.close('all')
        
    
    
    
        