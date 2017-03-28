execfile('go_lite')
import taxi
from p19b_particle_tools import *
#Clump pre-image analysis
#1.) Run amazing huge sims.  The best.
#2.) Run clumps on high density regions.  
#3.) Save Leaf Clumps.  (Probably with Britton's tool. )
# .) For each clump
#4.)   Find clump indices from peaks
#5.)   Find data on prior snapshots
#6.)     -- Projections for simple vis
#7.)     -- Fake herschel, alma maps; local (region) imges, full cloud images.
#8.)     -- Physical Properties.

def clump_finder(ds, loc = None, width = None):
    if loc is None:
        value, loc = ds.find_max('density')
    if width is None:
        width = (0.05,'code_length')
    sphere = ds.sphere(loc,width)
    master_clump = Clump(sphere,"density")
    master_clump.add_validator("min_cells", 20)
    #master_clump.add_validator("gravitationally_bound", use_particles=False, use_thermal_energy=False)
    c_min = sphere["gas", "density"].min()
    c_max = sphere["gas", "density"].max()
    step = 100
    find_clumps(master_clump, c_min, c_max, step)
    return sphere, master_clump, loc

if 'u05' not in dir():
    u05 = taxi.taxi('u05')
    car = u05

#
# Get clumps.
#
if 0:
    car.fill(125)
    ds = car.ds
    sphere, master_clump, loc = clump_finder(ds,width=(0.05,'code_length'), loc = ds.arr([ 0.03613281,  0.79589844,  0.03027344], 'code_length'))
    proj = ds.proj('density',2,data_source=sphere,center=loc)
    leaf_clumps = get_lowest_clumps(master_clump) #if both min_cells and grav_bound are used, this is empty.
    pw = proj.to_pw(center = loc, width = (0.1,'code_length'))
    pw.annotate_clumps(clump_list_sort(leaf_clumps))
    pw.save('clump_testx')
elif 1:
    car.fill(125)
    ds = car.ds
    value, loc = ds.find_max('density')
    min_dx = ds.index.get_smallest_dx()
    region = ds.region(center=loc,left_edge = loc-1.5*min_dx, right_edge = loc+1.5*min_dx)
    indices = region['particle_index'].astype('int64')
    leaf_indices=[indices]

#
# pre-images
#
if 0:
    for frame in [10]:
        car.fill(frame)
        ds  = car.ds
        proj_full = ds.proj('density',0,center='c')
        pw_full = proj_full.to_pw(center = 'c',width=(1.0,'code_length'))
        pw_full.set_cmap('density','gray')

        for nc,indices in enumerate(leaf_indices):
            pw_full.annotate_select_particles(1.0, col='r', indices=indices)
        pw_full.save('u05_peak_%04d'%frame)

if 1:
    for frame in [10]:
        car.fill(frame)
        ds  = car.ds
        ad = ds.all_data()
        for nc,indices in enumerate(leaf_indices):
            mask_to_get = np.zeros(indices.shape, dtype='int32') #this is necessary.  Factor of 2 faster to store mask_to_get
            t0=time.time()
            ad.set_field_parameter('indices_late',indices)
            ad.set_field_parameter('mask_to_get',mask_to_get)
            ad.set_field_parameter('timer',[t0])
            deposit_field = ad[("deposit",'deposit_target_particles')] #I believe this is necessary to build the field initially
            cut_region  = ad.cut_region(['obj["deposit","deposit_target_particles"] > 0'])
            mass = ( cut_region['density']*cut_region['cell_volume']).sum()
            print mass


