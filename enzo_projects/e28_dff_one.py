
dir1 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/MHDCTOrszagTangAMR'
dir2 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/MHDCTOrszagTangAMR-CenteredB'
#dir1 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/MHDCTOrszagTang'
#dir2 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/MHDCTOrszagTang-CenteredB'
#dir1 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/LoopAdvection_CT'
#dir2 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/LoopAdvection_CT_CenteredB'
##dir1 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/LoopAdvection_Dedner'
#dir2 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/LoopAdvection_Dedner_CenteredB'
#dir1 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/MHDCTOrszagTang-Restart'
#dir2 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/MHDCTOrszagTang-CenteredB'


cpu=0 #can be made fancier.
if 'extra' not in dir():
    frame = 0
    extra = -1
    grid1=1
    grid2=grid1
    field1 = 'Bx'
    field2 = field1
field2 = field1
grid2 = grid1
theset = 'DD%04d/data%04d'%(frame,frame)
if extra > 0:
    theset = 'ED%02d_%04d/Extra%02d_%04d'%(extra,frame,extra,frame)
fname1 = "%s/%s"%(dir1,theset)
fname2 = "%s/%s"%(dir2,theset)
cpuname1 = "%s.cpu%04d"%(fname1,cpu)
cpuname2 = "%s.cpu%04d"%(fname2,cpu)
print cpuname1
print cpuname2
fptr1 = h5py.File(cpuname1)
fptr2 = h5py.File(cpuname2)
try:
    grid_set1 = fptr1['Grid%08d'%grid1]
    grid_set2 = fptr2['Grid%08d'%grid2]
    field_set1 = grid_set1[field1][:]
    field_set2 = grid_set2[field2][:]

    if (field_set1.shape != field_set2.shape):
        print "set shape mismatch", field_set1.shape , field_set2.shape
    else:
        diff = field_set1 - field_set2
        stat(field_set1,'%s 1'%field1)
        stat(field_set2,'%s 2'%field2)
        stat(diff,'diff')
except:
    raise
finally:
    fptr1.close()
    fptr2.close()

