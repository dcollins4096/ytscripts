

dir1 = '/scratch1/dcollins/TestRunner/dev/5d6653715fb6_E28_enzo-dev/MHD/2D/MHDCTOrszagTang'
dir2 = '/scratch1/dcollins/TestRunner/dev/eae35fd882be+_E28_cbe_single_test_w04/MHD/2D/MHDCTOrszagTang'
frame = 4

ds1 = yt.load("%s/DD%04d/data%04d"%(dir1,frame,frame))
ds2 = yt.load("%s/DD%04d/data%04d"%(dir2,frame,frame))
field = 'density'
proj1 = yt.ProjectionPlot(ds1,2,field)
print proj1.save('OT1_%04d'%frame)
proj2 = yt.ProjectionPlot(ds2,2,field)
print proj2.save('OT2_%04d'%frame)
def gd(field):
    g1 = ds1.index.grids[0][field].v
    g2 = ds2.index.grids[0][field].v
    stat(g1-g2/(0.5*(g1+g2)))
