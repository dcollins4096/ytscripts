
if 'ef' not in dir():
    execfile('go')
import access_thing
reload(access_thing)
import clfind2d
import pyfits
def fake_temp(field,data):
    return data.ds.arr(np.ones_like(data['density'].v)*10,'K')
ef('p14_clump_mask.py')

yt.add_field('temperature',function=fake_temp,units='K')

fits_basedir='/scratch1/dcollins/Paper14b/'
at30 = access_thing.access_thing('high_256',30,fits_basedir=fits_basedir,frb_resolution=256)
at50 = access_thing.access_thing('high_256',50,fits_basedir=fits_basedir,frb_resolution=256)
at80 = access_thing.access_thing('high_256',80,fits_basedir=fits_basedir,frb_resolution=256)
at50L = access_thing.access_thing('low_256',50,fits_basedir=fits_basedir,frb_resolution=256)
at50L512 = access_thing.access_thing('low_512',50,fits_basedir=fits_basedir,frb_resolution=256)
at = at50L512
ef('unit_tool.py')
AvBar = L_cgs*n0_cgs * AvPerColumn
ThisAv = 3

ds = yt.load(at.enzo_dataset)
proj = yt.ProjectionPlot(ds,1,'density')
proj.set_cmap('density','gray')
proj.set_zlim('density',0.1/AvBar,50/AvBar)
proj.save('p14b_cmap1')
