import access_thing
import clfind2d
import pyfits
def fake_temp(field,data):
    return data.ds.arr(np.ones_like(data['density'].v)*10,'K')
yt.add_field('temperature',function=fake_temp,units='K')

fits_basedir='/scratch1/dcollins/Paper14b/'
if 'at' not in dir():
    reload(access_thing)
    at30 = access_thing.access_thing('high_256',30,fits_basedir=fits_basedir,frb_resolution=256)
    at50 = access_thing.access_thing('high_256',50,fits_basedir=fits_basedir,frb_resolution=256)
    at80 = access_thing.access_thing('high_256',80,fits_basedir=fits_basedir,frb_resolution=256)
    at50L = access_thing.access_thing('low_256',50,fits_basedir=fits_basedir,frb_resolution=256)
    at50L512 = access_thing.access_thing('low_512',50,fits_basedir=fits_basedir,frb_resolution=512)
    #at = at50
    at = at50L
    #at = at50L512
    #clump_name = 'cl2d_three_yt2'


if 1:
    todays_axis='x'
    ppv_name = 'ppv_4'
    ppv_args={'dims':256,'velocity_bounds':(-46,30,128,'code_velocity')}
    cube_read,cube_ds,vel = at.get_ppv(todays_axis,prefix=ppv_name, ppv_args=ppv_args)
