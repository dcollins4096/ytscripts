
#either with or without AMR
dir1 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/MHDCTOrszagTangAMR'
dir2 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/MHDCTOrszagTang-CenteredB'
dir1 = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/2D/MHDCTOrszagTang-Restart'

#frame = 2
#extra = 1
theset = 'DD%04d/data%04d'%(frame,frame)
#theset = 'ED%02d_%04d/Extra%02d_%04d'%(extra,frame,extra,frame)
print theset
ds1 = yt.load('%s/%s'%(dir1,theset))
ds2 = yt.load('%s/%s'%(dir2,theset))
proj1=yt.ProjectionPlot(ds1,2,'density')
proj1.annotate_grids()
proj1.save('OTAMRfid')
proj2=yt.ProjectionPlot(ds2,2,'density')
proj2.annotate_grids()
proj2.save('OTAMR_removal')
ng1= len(ds1.index.grids)
ng2 =len(ds2.index.grids)
if ng1 != ng2:
    print "ERROR different number of grids!", ng1, ng2
for g in range(ng1):
    ggg1 = ds1.index.grids[g]['density']
    ggg2 = ds2.index.grids[g]['density']
    print "total difference", g, np.abs(ggg1-ggg2).sum()
