
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
at = at50
#at.make_frb('density','x')
#d=at.get_dendrogram('x')
#density = at.get_fits_array('density','x')
#mask=at.get_cl2d('x')
#plt.clf()
#img=plt.imshow(np.log10(density),cmap='gray')
#plt.contour(mask)
#plt.savefig('p14b_contour1.png')
#at.image_cl2d( 'x','cl2d_two','p14b_cl2d_two_image')
#mask=at.get_cl2d('x','cl2d_two',clfind_args={'levels':[2],'log':True,'nPixMin':3})
#cube_read,cube_ds = at.get_ppv('x', 'ppv_1')
#masses_code = at.mass_list('x','cl2d_five',background_column=5)
#cube_read,cube_ds,vel = at.get_ppv('x',prefix='ppv_2', ppv_args={'dims':256,'velocity_bounds':(-35,35,128,'code_velocity')})
#at.spectra( 'x','cl2d_two','ppv_2')
#vel, spectra = at.spectra_image_pair( 'x','cl2d_two','ppv_1', image_zoom=False)
#at.image_cl2d( 'x','cl2d_two','ppv_1','p14b_cl2d_two_image')
#alpha = at.alpha( 'x','cl2d_two','ppv_2',background_column = 2)

if 0:
    ef('unit_tool.py')
    AvBar = L_cgs*n0_cgs * AvPerColumn
    ThisAv = 3
    at = at50L512,
    at.get_cl2d('x','cl2d_Av3',clfind_args={'levels':[ThisAv/AvBar,300], 'log':True, 'nPixMin':3})
    at.image_cl2d('x',clump_prefix='cl2d_Av3',out_prefix='p14b_Av3_cl2d')

if 0:

    at = at50L512
    at.pdf('x')
if 1:
    #mask=at.get_cl2d('x','cl2d_two')

    clfind_args = {'levels':[3], 'log':True,'nPixMin':3}
    clump_name = 'cl2d_three_yt3'
    mask=at.get_cl2d('x',clump_name,clfind_args=clfind_args)
    cube_read,cube_ds,vel = at.get_ppv('x',prefix='ppv_1')
    alpha = at.alpha( 'x',clump_name,'ppv_1',background_column = 3)
    masses = at.masses('x',clump_name, background_column = 3)
    #at.get_cut_region('x','cl2d_three_yt2', 3)
    at.get_clumps('x',clump_name,3) #clump = 3
    #at.compute_clump_properties()

#at.alpha_mass('x','cl2d_two','ppv_1',background_column=2, with_numbers=True)
if 0:
    #mask=at.get_cl2d('x','cl2d_two')
    mask=at.get_cl2d('x','cl2d_two')
    for clump in np.unique(mask)[1:]:
        at.get_clumps('x','cl2d_two',clump)
        at.compute_clump_properties()
    #at.get_clumps('x','cl2d_two',7)
    #at.compute_clump_properties()

if 0:
    ef('p14b_alpha.py')

#ef('p14b_chanels.py')
#cube_read,cube_ds,vel = at.get_ppv('x',prefix='ppv_1', ppv_args={'dims':256,'velocity_bounds':(-35,35,14,'code_velocity')})
#print np.isnan(cube_read).sum()
#at.image_cl2d( 'x','cl2d_two','p14b_cl2d_two_image')
