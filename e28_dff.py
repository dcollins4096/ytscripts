import h5py
dir1 = '/scratch1/dcollins/TestRunner/dev/0f26d0afc890+_E28_cbe_single_test_finally/MHD/2D/SedovBlast-MHD-2D-Fryxell'
dir2 = '/scratch1/dcollins/TestRunner/dev/cb9712a1f763+_enzo-dev_gold_push_with_no_specific/MHD/2D/SedovBlast-MHD-2D-Fryxell'
dir1 = '/scratch1/dcollins/TestRunner/dev/0f26d0afc890+_E28_cbe_single_test_finally/Hydro/Hydro-2D/RadiatingShockWave'
dir2 = '/scratch1/dcollins/TestRunner/dev/cb9712a1f763+_enzo-dev_gold_push_with_no_specific/Hydro/Hydro-2D/RadiatingShockWave'
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
    field1 = 'Density'
    field2 = field1
field2 = field1
grid2 = grid1
theset = 'DD%04d/data%04d'%(frame,frame)
theset = 'DD%04d/DD%04d'%(frame,frame)
if extra > 0:
    theset = 'ED%02d_%04d/Extra%02d_%04d'%(extra,frame,extra,frame)
fname1 = "%s/%s"%(dir1,theset)
fname2 = "%s/%s"%(dir2,theset)
if 1:
    ds1 = yt.load(fname1)
    ds2 = yt.load(fname2)
    g1 = ds1.index.grids[0]
    g2 = ds2.index.grids[0]
    field_set1 = g1[field1].v 
    field_set2 = g2[field2].v
    diff = (field_set1 - field_set2)
    fieldavg = (0.5*(field_set1 + field_set2))
    nonzero = diff != 0
    reldiff = diff[nonzero]/fieldavg[nonzero]
    
    stat(field_set1,'%s 1'%field1)
    stat(field_set2,'%s 2'%field2)
    stat(diff,'diff')
    stat(reldiff,'reldiff')

if 0:
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
            diff = (field_set1 - field_set2)
            fieldavg = (0.5*(field_set1 + field_set2))
            nonzero = diff != 0
            reldiff = diff[nonzero]/fieldavg[nonzero]
            
            stat(field_set1,'%s 1'%field1)
            stat(field_set2,'%s 2'%field2)
            stat(diff,'diff')
            stat(reldiff,'reldiff')
    except:
        raise
    finally:
        fptr1.close()
        fptr2.close()

