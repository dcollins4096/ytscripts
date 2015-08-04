

dir1 = '/scratch1/dcollins/EnzoProjects/E29_p10Boundary/ac01_unigrid_5step' 
dir2 = '/scratch1/dcollins/EnzoProjects/E29_p10Boundary/ac02_parallel_5step'
import dsdiff
reload(dsdiff)
frames=[0]
grids=[1]
ud = dsdiff.udiff(dir1,dir2,frames=frames,grids=grids, fields=['Density'])
##ud(output='y',raster='y',interogate_range=range(32))
ud(output='y', interogate_range=range(32))

