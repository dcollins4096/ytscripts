execfile('go')
from yt.analysis_modules.level_sets.api import * #for clumps
import pyximport; pyximport.install()
import particle_ops
import clump_particles
reload(clump_particles)
thisdir = '/scratch1/dcollins/Paper19/B02/u02-128-d-alpha0.5'
template = thisdir+'/DD%04d/data%04d'
ef('particle_selector.py')

if 0:
    """Make a sphere"""
    thisdir = '/scratch1/dcollins/Paper19/B02/u02-128-d-alpha0.5'
    frame = 360
    width = (0.05,'code_length')
    setname = template%(frame,frame)
    ds = yt.load(setname)
    width = ds.arr(0.05,'code_length')
    val,loc = ds.find_max('density')
    #loc = ds.arr([0.5]*3,'code_length')
    """Image that sphere"""
    sphere = ds.sphere(loc,width)
    proj2 = ds.proj('density',0,data_source=sphere,center=loc)
    pw = proj2.to_pw(center = loc, width = (0.1,'code_length'))
    pw.save('test2.png')

if 0:
    """Full image, just to check."""
    proj_full = ds.proj('density',0, center = loc) #width = (1.0,'code_length'))
    pw_full = proj_full.to_pw(center = loc,width=(1.0,'code_length'))
    pw_full.annotate_particles(1.0,stride=16)
    pw_full.save('test_full.png')

if 0:
    """Get the indices.  Save them with a pickle."""
    indices_late,xpos_late,ypos_late,zpos_late = clump_particles.particles_from_clump(sphere)
    #fPickle.dump(indices_late,'p19_u02d_0360_indices_late.pickle')

if 0:
    """Double check"""
    pw.annotate_dave_particles(1.0, indices=indices_late)
    pw.save('test3.png')
   
if 0:
    """From the loaded particles, make a test sphere"""
    indices_0360 = fPickle.load('p19_u02d_0360_indices_late.pickle')
    frame = 360
    setname = template%(frame,frame)
    ds_0360 = yt.load(setname)
    loc = ds_0360.arr([ 0.03613281,  0.79589844,  0.03027344], 'code_length') #the max
    width = ds_0360.arr(0.05,'code_length')
    sphere = ds_0360.sphere(loc,width)


if 0:
    """Get deposit_target_particles_1;  both field parameters are necessary."""
    t0=time.time()
    print "sum indices 0360", indices_0360.sum()
    #indices_0360.sort() #this is ok, this list doesn't have any order dependancies.  Not necssary now, later maybe.
    mask_to_get = np.zeros(indices_0360.shape, dtype='int32') #this is necessary.  Factor of 2 faster to store mask_to_get
    sphere.set_field_parameter('indices_late',indices_0360)
    sphere.set_field_parameter('mask_to_get',mask_to_get)
    tdep = sphere[("deposit",'deposit_target_particles_1')]
    dep = sphere[("deposit",'all_count')]
    t1 = time.time()
    print "zero is ok", np.abs(tdep-dep).sum(), "dt = ", t1-t0

if 0:
    """ make a movie with the particles.  Stride 64 over particles to make it go."""
    indices_0360 = fPickle.load('p19_u02d_0360_indices_late.pickle')
    for frame in [360]:# [340,350, 360]:
        setname = template%(frame,frame)
        ds = yt.load(setname)
        proj_full = ds.proj('density',0, center = 'c' ) #width = (1.0,'code_length'))
        pw_full = proj_full.to_pw(center = 'c',width=(1.0,'code_length'))
        pw_full.annotate_dave_particles(1.0, indices=indices_0360[::64])
        pw_full.save('u02d_history_%04d'%frame)

if 0:
    """Movie of the deposited field."""
    indices_0360 = fPickle.load('p19_u02d_0360_indices_late.pickle')
    indices_0360 = indices_0360
    t0 = time.time()
    for frame in [360]:# [340,350, 360]:
        setname = template%(frame,frame)
        ds = yt.load(setname)
        ad = ds.all_data()
        mask_to_get = np.zeros(indices_0360.shape, dtype='int32') #this is necessary.  Factor of 2 faster to store mask_to_get
        ad.set_field_parameter('indices_late',indices_0360)
        ad.set_field_parameter('mask_to_get',mask_to_get)
        ad.set_field_parameter('timer',[t0])
        proj_full = ds.proj(("deposit",'deposit_target_particles_1'),0, data_source=ad, center = 'c' ) #width = (1.0,'code_length'))
        pw_full = proj_full.to_pw(center = 'c',width=(1.0,'code_length'))
        pw_full.save('u02d_deposited_%04d'%frame)
    t1 = time.time()
    print "time = ", t1-t0

if 0:
    frame = 360
    setname = template%(frame,frame)
    """ now let's try out the clump tools"""
    frame = 360
    setname = template%(frame,frame)
    ds = yt.load(setname)
    #val,loc = ds.find_max('density') (the max, because I'm super lazy)
    loc = ds.arr([ 0.03613281,  0.79589844,  0.03027344], 'code_length')
    width = (0.05,'code_length')
    sphere = ds.sphere(loc,width)
    master_clump = Clump(sphere,"density")
    master_clump.add_validator("min_cells", 20)
    #master_clump.add_validator("gravitationally_bound", use_particles=False, use_thermal_energy=False)
    c_min = sphere["gas", "density"].min()
    c_max = sphere["gas", "density"].max()
    step = 10
    find_clumps(master_clump, c_min, c_max, step)
    leaf_clumps = get_lowest_clumps(master_clump) #if both min_cells and grav_bound are used, this is empty.
    proj2 = ds.proj('density',2,data_source=sphere,center=loc)
    pw = proj2.to_pw(center = loc, width = (0.1,'code_length'))
    pw.annotate_clumps(leaf_clumps)
    pw.save('clump_testx')

if 0:
    for n, c in enumerate(leaf_clumps):
        print c['density'].shape, c['cell_mass'].shape, c.quantities.total_quantity('cell_mass') #, c.quantities.total_mass()

if 0:
    import clump_stuff
    reload(clump_stuff)
    #for c in leaf_clumps

