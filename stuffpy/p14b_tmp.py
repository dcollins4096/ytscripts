
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
at = access_thing.access_thing('high_256',30,fits_basedir=fits_basedir,frb_resolution=256)
clump_prefix='cl2d_two'
axis='x'
twod_clump = np.transpose(at.get_cl2d(axis,clump_prefix))
flat_clump = twod_clump.flatten()
clump_stuff_tuple = (flat_clump,'xyz'.index(axis),twod_clump)
clump_stuff_dict={'clump_mask_stuff':clump_stuff_tuple}
ds = yt.load(at.enzo_dataset)
ef('p14_get_cut_region.py')

cr = at.get_cut_region(axis, clump_prefix, clump_ind)


