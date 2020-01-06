from go import *


car = taxi.load('zc02')
car.fields=['density','x-velocity','y-velocity','z-velocity','magnetic_field_x',\
            'magnetic_field_y','magnetic_field_z']
car.frames=[100]
car.axis=[0,1,2]
#car.plot()
if 0:
    for frame in car.frames:
        ds = car.load(frame)
        level=0
        dims=[256]*3
        field_list=car.fields
        cg = ds.covering_grid(level,left_edge=[0.0,0.0,0.0],dims=dims,fields=field_list)
        for field in field_list:
            outfile = "p57_cubes/cube_%s_%04d_%s.fits"%(car.outname,frame,field)
            hdu = pyfits.PrimaryHDU(cg[field])
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto(outfile,clobber=True)

if 1:
    for frame in car.frames:
        ds = car.load(frame)
        level=0
        dims=[256]*3
        field_list=car.fields
        cg = ds.covering_grid(level,left_edge=[0.0,0.0,0.0],dims=dims,fields=field_list)
        for field in field_list:
            outfile = "p57_cubes/cube_%s_%04d_%s.fits"%(car.outname,frame,field)
            hdu = pyfits.PrimaryHDU(cg[field])
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto(outfile,clobber=True)

