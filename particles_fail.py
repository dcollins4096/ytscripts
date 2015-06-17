import yt
from yt.analysis_modules.level_sets.api import  *
import numpy as np

if 0:
    setname = '/scratch1/dcollins/Paper19/ClumpTest/IsolatedGalaxy_Gravity/galaxy0030/galaxy0030'
    ds1 = yt.load(setname)
    step1 = 10 

    """this doesn't work."""
    ds2 = yt.testing.fake_amr_ds( fields = ("particle_position_x",
                                         "particle_position_y",
                                         "particle_position_z",
                                         "particle_mass",
                                         "particle_velocity_x",
                                         "particle_velocity_y",
                                         "particle_velocity_z", "density"))
    step2 = 1+5e-1

    #ds1, step1 fails.  ds2, step2 does make a bunch of clumps, but doesn't seem to fail.  
    ds = ds1
    step = step1
    val,loc = ds.find_max('density') 
    width = (0.05,'code_length')
    sphere = ds.sphere(loc,width)
    master_clump = Clump(sphere,"density")
    master_clump.add_validator("min_cells", 20)
    c_min = sphere["gas", "density"].min()
    c_max = sphere["gas", "density"].max()
    n_steps = (np.log(c_max)-np.log(c_min))/np.log(step)
    print n_steps

    find_clumps(master_clump, c_min, c_max, step)
    leaf_clumps = get_lowest_clumps(master_clump) #if both min_cells and grav_bound are used, this is empty.
    print "leaf clumps", leaf_clumps

    for n, c in enumerate(leaf_clumps):
        print c['density'].shape, c['cell_mass'].shape, c.quantities.total_quantity('cell_mass'), c.quantities.total_mass()
