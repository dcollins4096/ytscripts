if 'ef' not in dir():
    execfile('go')
def code_mag_eng(field,data):
    bx = data['Bx'].in_units('code_magnetic')
    by = data['By'].in_units('code_magnetic')
    bz = data['Bz'].in_units('code_magnetic')
    out = 0.5*(bx*bx+by*by+bz*bz).v
    return out
yt.add_field('code_mag_eng',code_mag_eng)
def code_db_over_b(field,data):
    return data['DivB'].v/data['code_mag_eng']
yt.add_field('code_db_over_b',code_db_over_b)

basedir = '/Users/dcollins/scratch/EnzoProjects/E12/'

c01 = 'c01_ppm'
c02 = 'c02_ct'
c03 = 'c03_ct_bsmall'
c04 = 'c04_ppm_drill'
c05 = 'c05_ct_small_betterCT'
c06 = 'c06_c05_better_divb'
c08 = 'c08_dev_ctBalsara'
c09 = 'c09_dev_ct_details' #This one is a full cosmology run.  Results in the run dir/images
r01 = 'r01_turb_newct'
s01 = 's01_dev_sphere1'
set1 = s01
dirname = "%s/%s"%(basedir,set1)
if 'frame' not in dir():
    frame = 2
if 'divb' not in dir():
    divb = []
if 'scales' not in dir():
    scales = []
ds_name = dsdiff_helpers.get_ds_name(dirname,frame)
#ds_name = "%s/DD%04d/DD%04d"%(dirname,frame,frame)

ds = yt.load(ds_name)
if 0:
    z = ds['CosmologyCurrentRedshift']
    a = 1/(1+z)
    scales.append([frame,a])
ad = ds.all_data()
stat(ad['DivB'], "divb frame %d"%frame)
#divb.append([frame,np.abs(ad['DivB']).max()])
if 'axis' not in dir():
    axis = 0

fields = ['grid_level', 'DivB','density', 'magnetic_energy']
#fields = ['density', 'grid_level']
#fields = []
#fields = ['DivB']
#fields = ['grid_level']
#fields=[]
for field in fields:
    proj = yt.ProjectionPlot(ds,axis,field)
    proj.annotate_grids()
    print proj.save('%s_%04d'%(set1,frame))
#ds.print_stats()
