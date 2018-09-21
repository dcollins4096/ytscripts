

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
w03 = 'w03_dev_amr'
w04 = 'w04_cbe_amr'
set1 = w03
set2 = w04
dir1 = "%s/%s"%(basedir,set1)
dir2 = "%s/%s"%(basedir,set2)
dir1 = '/scratch1/dcollins/TestRunner/dev/0f26d0afc890+_E28_cbe_single_test_finally/Hydro/Hydro-2D/RadiatingShockWave'
dir2 = '/scratch1/dcollins/TestRunner/dev/cb9712a1f763+_enzo-dev_gold_push_with_no_specific/Hydro/Hydro-2D/RadiatingShockWave'
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

