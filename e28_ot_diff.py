

ot3 = 'ot3_dev'
ot1 = 'o1_CB'
#ot2 = '/Users/dcollins/scratch/TestRunner/920fe8acac49+_CenteredBremoval_arse/MHD/2D/MHDCTOrszagTang'
ot2 = 'o2_dve'
ot4 = 'ot4_1dreprise'
basedir = '/Users/dcollins/scratch/EnzoProjects/E28_CenteredBRemoval/OT_b/'
ot4 = 'ot4_1dreprise'
ota = 'ota_badone'
otb = 'otb_failure'
otc = 'otc_otb2'
w01 = 'w01_enzo_dev'
w02 = 'w02_CBE'
set1 = w01
set2 = w02
dir1 = "%s/%s"%(basedir,set1)
dir2 = "%s/%s"%(basedir,set2)
import dsdiff
reload(dsdiff)
frames=[0]
grids=[[1,1]] #,[1,2]]

print "\n\n\n"
#ds1 = yt.load(dir1+"/DD%04d/data%04d"%(frames,frames))
#ds2 = yt.load(dir2+"/DD%04d/data%04d"%(frames,frames))
#field = 'Bx'
#print yt.ProjectionPlot(ds1,2,field).save(set1)
#print yt.ProjectionPlot(ds2,2,field).save(set2)
ud = dsdiff.udiff(dir1,dir2,frames=frames,grids=grids, fields=['Density'])
ud()
##ud(output='y',raster='y',interogate_range=range(32))
#ud(output='y', interogate_range=range(32))
