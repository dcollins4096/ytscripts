import pyfits
base = '/scratch1/dcollins/Paper08/B02/512'
if 'frame' not in dir():
    frame = 0
pfname = '%s/RS%04d/restart%04d'%(base,frame,frame)

outname = 'High_512'
size = [512]*2
ax = 0
pf = yt.load(pfname)
cg = pf.covering_grid(0,[0.0]*3, [512]*3, fields=['x-velocity','y-velocity','z-velocity','Density'])

#proj = pf.proj('Density',ax)
#frb = proj.to_frb((1.0, 'code_length'),size)
for field in ['x-velocity','y-velocity','z-velocity','Density']:
    setname = "%s_%04d_%s_cube.fits"%(outname,frame,field)
    #thisset = frb[field]
    hdu = pyfits.PrimaryHDU(cg[field])
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(setname)
    hdulist.close()
    print setname
