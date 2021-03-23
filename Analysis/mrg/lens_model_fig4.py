
def make_fig4():
    
    from grizli.pipeline import auto_script
    from mastquery import overlaps
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
    import numpy as np      
    
    from grizli import utils
    
    root = 'j013804m2156'

    #x0, y0= 2465, 1740
    x0, y0 = 1906, 1358
    
    NY = (2204-y0)//2
    NX = (2898-x0)//2
    NX = 520
    
    x0 += NX
    y0 += NY
    
    slx = slice(x0-NX, x0+NX)
    sly = slice(y0-NY, y0+NY)
    
    mi = -0.01
    
    si = 4.8
    si = 8
    
    dpi = 100
    
    for i in range(3):
        _rgb = auto_script.field_rgb(root=root, filters=['f105w','f125w','f160w'], HOME_PATH=None, xyslice=(slx, sly), show_ir=False, add_labels=False, output_format='png', suffix='.lensmodel.on', rgb_min=mi, tick_interval=10, xsize=si, output_dpi=dpi) #int(2*NX/si))
        
        fig = _rgb[-1]
        ax = fig.axes[0]
        #ax.axis('off')
        fig.tight_layout(pad=0)
        fig.savefig('j013804m2156.lensmodel.on.png')
        
        img = plt.imread('j013804m2156.lensmodel.on.png')
        sh = img.shape
        dpi = int(np.round(dpi*sh[0]/2/NY))
        print(i, dpi, sh, 2*NX, 2*NY)
        
    plt.close('all')
    
    import pyregion
    #r2 = pyregion.parse(region_string).as_imagecoord(f[0].header)
    #patch_list, artist_list = r2.get_mpl_patches_texts()
    
    im = pyfits.open(f'{root}-f160w_drz_sci.fits')
    wcs = pywcs.WCS(im[0].header, relax=True)
    wsl = wcs.slice((sly, slx))
    hsl = utils.get_wcs_slice_header(wcs, slx, sly)
    
    reg = pyregion.open('mrg0138_supernova/Analysis/lensing/johan_allpot_adjust.reg').as_imagecoord(hsl)
    #reg_sky = pyregion.open('johan_allpot.reg').as_imagecoord(hsl)
    
    img = plt.imread('j013804m2156.lensmodel.on.png')
    sh = img.shape
    
    img = img[::-1,:,:]
    
    sx = 8

    patches, artists = reg.get_mpl_patches_texts()

    fig, ax = plt.subplots(1,1,figsize=[sx, sx*sh[0]/sh[1]])
    ax.imshow(img, origin='lower') #, extent=wsl.calc_footprint()[[0,2],:].T.flatten())
    
    for pat in patches:
        pat.set_zorder(100)
        ax.add_patch(pat)
        
    for txt in artists:
        txt.set_zorder(100)
        if 0:
            ax.add_artist(txt)
    
    ########## Labels
    coo = utils.read_catalog('coords.txt')
    N, inter = 10, 0.5
    N, inter = 20, 1.0
    xy0 = np.array(wcs.all_world2pix(coo['ra'], coo['dec'], 0)).T
    xy = np.cast[int](np.round(wcs.all_world2pix(coo['ra'], coo['dec'], 0))).T
    xyi = (xy - np.array([x0-NX, y0-NY]))/np.array([NX, NY])/2*np.array(sh[:-1])
    xyi0 = (xy0 - np.array([x0-NX, y0-NY]))/np.array([NX, NY])/2*np.array(sh[:-1])
    #dy = 20
    pixgrow = sh[0]/2/NY
    
    xyi = np.array(wsl.all_world2pix(coo['ra'], coo['dec'], 0)).T
    
    # Label offsets
    dx = [-25, 2, 0, -28, 0, 0, -28, -27] + [22, 22, -22, 3]
    dy = [-8, 18, 12, 0, 10, 10, -15, 10] + [2, 2, 0, 10]
    
    img0 = plt.imread('j013804m2156.on.png')
    sh0 = img.shape
    dx = np.array(dx)*sh[1]/sh0[1]
    dy = np.array(dy)*sh[0]/sh0[0]
    
    colors = ['w']*7 + ['pink'] + ['w']*4

    for j, i in enumerate(range(8+4)):#[4,5,6]):            
        # if i == 7:
        #     col = 'pink'
        # elif i > 7:
        #     col = '0.5'
        # else:
        #     col = 'w'
        
        if i == 7:
            continue
            
        if i > 7:
            fs = 8
        else:
            fs = 8
        
        ax.scatter(xyi[i][0], xyi[i][1], marker='x', color=colors[i])
            
        ax.text(xyi[i][0]+dx[i], xyi[i][1]-dy[i], coo['label'][i], ha='center', va='top', color=colors[i], fontsize=fs)
    
    
    ax.axis('off')
    fig.tight_layout(pad=0)
    
    fig.savefig('fig4_lensmodel_layout.pdf')
    