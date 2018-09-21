

ac01n0 = '/scratch1/dcollins/EnzoProjects/E29_p10Boundary/ac01_unigrid_5step' 
ac02n0 = '/scratch1/dcollins/EnzoProjects/E29_p10Boundary/ac02_parallel_5step'
ac01='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac01_uni'
ac02 = '/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac02_2core'

ac03n='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/FromNazare/ac03_ac01again'
ac01n='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/FromNazare/ac01_unigrid_5step'
ac04n='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/FromNazare/ac04_ac02again'
ac02n='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/FromNazare/ac02_parallel_5step'
ac07 ='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac07_ot_serial'
ac08 ='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac08_ot_parallel'
ac09 ='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac09_uni_periodic'
ac10 ='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac10_par_periodic'
ac11 ='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac11_serial_turb'
ac12 ='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac12_par_turb'
ac15 ='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac15_ser_simpleics'
ac16 ='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac16_par_simpleics'
ac17 ='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac17_ser_simple_32'
ac18 ='/Users/dcollins/scratch/EnzoProjects/E29_p10boundary/ac18_par_simple_32'
dir1 = ac15
dir2 = ac16
import dsdiff
reload(dsdiff)
frames=[5]
grids=[[1,1]] #,[1,2]]

print "\n\n\n"
ud = dsdiff.udiff(dir1,dir2,frames=frames,grids=grids, fields=['Density'])
ud(shift=[0.0, 0.0,0.0])
##ud(output='y',raster='y',interogate_range=range(32))
#ud(output='y', interogate_range=range(32))

