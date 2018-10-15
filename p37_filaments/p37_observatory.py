
ef('p37_blur.py')

if 1:
    if 'frame' not in dir():
        frame = 70

    ds = yt.load('/scratch1/dcollins/Paper08/B02/512/RS%04d/restart%04d'%(frame,frame)); simname = 'high_field'
    field = 'density'
    proj = ds.proj(field,"z")
    pw=proj.to_pw()

if 1:
    width = 1.0
    Nzones = 8192
    frb = proj.to_frb( width, int(Nzones))
    #plt.clf()
    #plt.imshow(np.log10(frb['density']), origin='lower',interpolation='nearest',cmap='gray')
    den = frb['density'].v
    distance = 200. #pc
    resolution = 18 #40.#arcsec
    smooth_pc = distance*resolution/206264. #in pc
    smooth_px = int(smooth_pc/(4.6/Nzones))
    blur = blur_image(den,smooth_px)
    hdu = pyfits.PrimaryHDU(blur)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('b02_512_smoothed_%04dzones_%fpc_%farcsec_density.fits'%(Nzones,distance,resolution))
    hdulist.close()
