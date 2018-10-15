#Clump pre-image analysis
#0.) Run amazing huge sims.  The best.
#1.) Run clumps on high density regions.  
#2.) Save Leaf Clumps.  (Probably with Britton's tool. )
# .) For each clump
#3.)   Find clump indices from peaks
#4.)   Find data on prior snapshots
#5.)     -- Projections for simple vis
#6.)     -- Fake herschel, alma maps; local (region) imges, full cloud images.
#7.)     -- Physical Properties.

def clump_finder(ds, loc = None, width = None):
    """step 1"""
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

if 'car' not in dir():
    u05 = taxi.taxi('u05')
    car = u05

#
# Get clumps.
#
if 'leaf_clumps' not in dir():
    """Step 3, version 2"""
    car.fill(frame)
    ds = car.ds
    sphere, master_clump, loc = clump_finder(ds,width=(0.05,'code_length'), loc = None ) #ds.arr([ 0.03613281,  0.79589844,  0.03027344], 'code_length'))
    proj = ds.proj('density',2,data_source=sphere,center=loc)
    leaf_clumps = get_lowest_clumps(master_clump) #if both min_cells and grav_bound are used, this is empty.
    pw = proj.to_pw(center = loc, width = (0.1,'code_length'))
    pw.annotate_clumps(clump_list_sort(leaf_clumps))
    pw.save('clump_testx')

#
# Get particle indices from said leaf clumps.
#
whole_clump_particles = False
if 'indices_late' not in dir() and whole_clump_particles:
    """this is broken"""
    from yt.analysis_modules.level_sets.api import * #for clumps
    import pyximport; pyximport.install()
    import particle_ops
    import clump_particles
    from yt.utilities.data_point_utilities import FindBindingEnergy
    reload(clump_particles)
    indices_late,xpos_late,ypos_late,zpos_late = clump_particles.particles_from_clump(sphere)

if 'indices_late' not in dir():
    peak_list = []
    for cl in leaf_clumps:
        peak_ind = np.where( cl['density'] == cl['density'].max())[0][0]
        xpeak = cl['x'][peak_ind]
        ypeak = cl['y'][peak_ind]
        zpeak = cl['z'][peak_ind]
        peak_list.append( nar([xpeak,ypeak,zpeak]))
    proj = yt.ProjectionPlot(car.ds,0,'density')
    for n,peak in enumerate(peak_list):
        proj.annotate_point(peak,n,text_args={'color':'r'})
    proj.set_cmap('density','Greys')
    proj.save('p19_%s_n%04d_peaks.png'%(car.outname,frame))
    fPickle.dump(peak_list,"p19_%s_512_n%04d_peaks2.pickle"%(car.outname,frame))
