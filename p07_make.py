import pyfits
base = '/scratch1/dcollins/Paper08/B02/512'
frame = 87
pfname = '%s/RS%04d/restart%04d'%(base,frame,frame)

outname = 'High_512'
size = [512]*2
ax = 0
pf = yt.load(pfname)
proj = pf.proj(field,ax)
frb = proj.to_frb((1.0, 'code_length'),size)
for field in ['By','Bz', 'Bx','Density']:
    setname = "%s_%04d_%s_%s.fits"%(outname,frame,field,'xyz'[ax])
    thisset = frb[field]
    hdu = pyfits.PrimaryHDU(thisset)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(setname)
    hdulist.close()
