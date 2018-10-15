
from yt.analysis_modules.level_sets.api import * #for clumps
dir1 = '/Users/dcollins/scratch/Paper05/OK4'
dir1 = '/scratch1/dcollins/Paper05/OK4/'
frame=700
size = 128*2**4
ds1 = yt.load(dir1+"/RS%04d/restart%04d"%(frame,frame))
ad = ds1.all_data()
master_clump = Clump(ad,"density")
master_clump.add_validator("min_cells", 20)
c_min = ad["gas", "density"].min()
c_max = ad["gas", "density"].max()
n_steps = 5
step_size = np.exp( (np.log(c_max) - np.log(c_min))/n_steps)
find_clumps(master_clump, c_min, c_max, step_size)
