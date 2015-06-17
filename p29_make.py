import pyfits
base = '/scratch1/dcollins/Paper08/B2/512'; outname = 'Mid_512'
base = '/scratch1/dcollins/Paper08/B02/512'; outname = 'High_512'
if 'frame' not in dir():
    frame = 2
pfname = '%s/RS%04d/restart%04d'%(base,frame,frame)

size1 = 512
ax = 0
pf = yt.load(pfname)
cg = pf.covering_grid(0,[0.0]*3, [size1]*3, fields=['x-velocity','y-velocity','z-velocity','Density'])

#proj = pf.proj('Density',ax)
#frb = proj.to_frb((1.0, 'code_length'),size)
for field in ['x-velocity','y-velocity','z-velocity','Density']:
    setname = "%s_%04d_%s_s_%d_cube.fits"%(outname,frame,field,size1)
    #thisset = frb[field]
    hdu = pyfits.PrimaryHDU(cg[field])
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(setname)
    hdulist.close()
    print setname
