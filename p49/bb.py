import astropy.io.fits as pyfits
frame=900
directory='/home/dcollins/scratch/p49/bb_01/256/b1p1/faun.rc.fas.harvard.edu/bburkhart/MHD/256/b1p1/t_%d'%frame
fname = {}
fname['density']=directory+"/dens_t%d.fits"%frame
fname['magnetic_field_x']=directory+"/magx_t%d.fits"%frame
fname['magnetic_field_y']=directory+"/magy_t%d.fits"%frame
fname['magnetic_field_z']=directory+"/magz_t%d.fits"%frame
fname['velocity_x']=directory+"/velx_t%d.fits"%frame
fname['velocity_y']=directory+"/vely_t%d.fits"%frame
fname['velocity_z']=directory+"/velz_t%d.fits"%frame
fname['Bx']=fname['magnetic_field_x']
fname['By']=fname['magnetic_field_y']
fname['Bz']=fname['magnetic_field_z']


data ={}
for f in fname:
    data[f]= np.array(pyfits.open(fname[f])[0].data)
bbox = np.array([[0.,1.]]*3) 
ds = yt.load_uniform_grid(data,data['density'].shape,length_unit='cm',bbox=bbox)
import turb_quan
reload(turb_quan)
qb = turb_quan.quan_box()
qb(ds=ds, frames=[frame])
