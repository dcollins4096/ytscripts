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
car.fill(125)
ds = car.ds
sphere, master_clump, loc = clump_finder(ds,width=(0.05,'code_length'), loc = ds.arr([ 0.03613281,  0.79589844,  0.03027344], 'code_length'))
proj = ds.proj('density',2,data_source=sphere,center=loc)
leaf_clumps = get_lowest_clumps(master_clump) #if both min_cells and grav_bound are used, this is empty.
pw = proj.to_pw(center = loc, width = (0.1,'code_length'))
pw.annotate_clumps(clump_list_sort(leaf_clumps))
pw.save('clump_testx')

#
# Get particle indices from said leaf clumps.
#
