
from yt.analysis_modules.level_sets.api import * #for clumps
frame = 125
fname = '/work/00369/tg456484/maverick/Paper19/u05-r4-l4-128/DD%04d/data%04d'%(frame,frame)

ds = yt.load(fname)

if 'loc' not in dir():
  val, loc = ds.find_max('density')

if 0:
  for ax in 'xyz':
    print "proj", ax
    proj = yt.ProjectionPlot(ds,ax,'density', center=loc)
    #proj.annotate_grids()
    print proj.save('u05_0125_nogrids')

if 0:
  prof = yt.ProfilePlot(ds.all_data(),'density','cell_volume',weight_field=None)
  prof.save('u05_0125')

if 1:
  master_clump = Clump(ds.all_data(),"density")
  master_clump.add_validator("min_cells", 20)
  c_min = 10
  c_max = 1e7
  step = 10
  find_clumps(master_clump, c_min, c_max, step)

