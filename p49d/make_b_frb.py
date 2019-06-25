import os
import sys
sys.path.append("%s/yt3_scripts"%os.environ['HOME'])
from go import *
import mpi4py
yt.enable_parallelism()
carname = sys.argv[1]
car = taxi.load(carname)
outputdir = "%s/frbs"%car.directory
clobber=False
my_axis = 'xyz'
#car.frames=[100]
#my_axis = car.axis
for frame in car.return_frames():
    ds = car.load(frame)
    for field in ['magnetic_field_strength']:
        for axis in my_axis:
            field_name = field + "_"+axis
            outfile = outputdir+"/DD%.4d_%s.fits" %(frame,field_name)
            if os.access(outfile, os.F_OK) and not clobber:
                print("FRB exists: %s"%outfile)
            else:
                proj = ds.proj(field,axis)
                res = ds.parameters['TopGridDimensions'][2 + ord('x') - ord(axis)] # zyx order
                frb = proj.to_frb(1,res)
                if yt.is_root():
                    hdu = pyfits.PrimaryHDU(frb[field])
                    hdulist = pyfits.HDUList([hdu])
                    hdulist.writeto(outfile,clobber=True)
                    print("wrote", outfile)

