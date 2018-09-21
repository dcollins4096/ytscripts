import numpy as na
import pdb
from dsdiff_helpers import *
class OverlapException(Exception):
    def __init__(self,n1,n2,grid1,grid2,field,left1,left2,right1,right2):
        self.value="No overlap: n(%d,%d) g(%d,%d), field %s"%(n1,n2,grid1,grid2,field)
    def __str__(self):
        return repr(self.value)
def  nint(x):
    return (x+0.5).astype('int')
    
def FindOverlap(dir1,dir2,n1,n2,grid1,grid2,field,shift=nar([0.0]*3),grid_direct=False,
               ds1=None,ds2=None, nGhostSkip=0):
    """Given two grids, return two slices that will index the overlap regions.
    *field* determines if overlap is considered for faces only.
    *shift* is applied to grid2, for periodic shifts
    OverlapException is raised if there's no overlap."""
    ds_name_1 = get_ds_name(dir1,n1)
    ds_name_2 = get_ds_name(dir2,n2)
    fake_grid_1 = fake_grid(dir1,n1,grid1)
    fake_grid_2 = fake_grid(dir2,n2,grid2)
    debug=-7

    
    left1  = fake_grid_1.GridLeftEdge+ shift
    right1 = fake_grid_1.GridRightEdge+ shift
    width1 = fake_grid_1.CellWidth
    left2  = fake_grid_2.GridLeftEdge
    right2 = fake_grid_2.GridRightEdge
    width2 = fake_grid_2.CellWidth
    #print "width1", width1
    #print "width2", width2

    if ds1 is not None:
        nGhost1 = ds1['NumberOfGhostZones']
        nGhost1 *= ds1['WriteBoundary']
        nGhost1 *= ds1['WriteGhostZones']
    if ds2 is not None:
        nGhost2 = ds2['NumberOfGhostZones']
        nGhost2 *= ds2['WriteBoundary']
        nGhost2 *= ds2['WriteGhostZones']
    if debug > 0:
        print "CLOWN left1", left1
        print "CLOWN left2", left2
        print "CLOWN right1", right1
        print "CLOWN right2", right2
    left1 -= nGhost1*width1; left2-=nGhost2*width2
    right1 += nGhost1*width1; right2+=nGhost2*width2
    if debug>0:
        print "CLOWN left1", left1
        print "CLOWN left2", left2
        print "CLOWN right1", right1
        print "CLOWN right2", right2

    strict_overlap_bool = [ left1 < right2 , left2 < right1 ]
    proper_overlap = strict_overlap_bool[0].all() and strict_overlap_bool[1].all()

    edge_overlap = [na.abs(left1- right2)< 0.1*width1, na.abs(left2 - right1) < 0.1*width2]
    #print "edge overlap", edge_overlap


    if not proper_overlap:
        #then there's a face overlap, whence check face/field pair

        #only the 6 listed fields can have an edge overlap
        if field not in  ['MagneticField_F_1','MagneticField_F_2','MagneticField_F_3']+\
               ['ElectricField_1','ElectricField_2','ElectricField_3']:
            raise OverlapException(n1,n2,grid1,grid2,field,left1,left2,right2,right2)  
        
        #x face overlap and Bx, Ey, Ez
        if (edge_overlap[0][0] or edge_overlap[1][0]) and \
            field not in ['MagneticField_F_1','ElectricField_2','ElectricField_3']:
            raise OverlapException(n1,n2,grid1,grid2,field,left1,left2,right2,right2)              
        #y face overlap and By, Ex, Ez
        if (edge_overlap[0][1] or edge_overlap[1][1]) and \
            field not in ['MagneticField_F_2','ElectricField_1','ElectricField_3']:
            raise OverlapException(n1,n2,grid1,grid2,field,left1,left2,right2,right2)              
        #z face overlap and Bz, Ex, Ey
        if (edge_overlap[0][2] or edge_overlap[1][2]) and \
            field not in ['MagneticField_F_3','ElectricField_1','ElectricField_2']:
            raise OverlapException(n1,n2,grid1,grid2,field,left1,left2,right2,right2)              
        
        #for more than one overlap, only the electric field.
        if (edge_overlap[0][0] or edge_overlap[1][0]) and \
           (edge_overlap[0][1] or edge_overlap[1][1]) and field != 'ElectricField_3':
            raise OverlapException(n1,n2,grid1,grid2,field,left1,left2,right2,right2)                          
        if (edge_overlap[0][0] or edge_overlap[1][0]) and \
           (edge_overlap[0][2] or edge_overlap[1][2]) and field != 'ElectricField_2':
            raise OverlapException(n1,n2,grid1,grid2,field,left1,left2,right2,right2)                          
        if (edge_overlap[0][2] or edge_overlap[1][2]) and \
           (edge_overlap[0][1] or edge_overlap[1][1]) and field != 'ElectricField_1':
            raise OverlapException(n1,n2,grid1,grid2,field,left1,left2,right2,right2)                          


    OverlapLeft =  na.amax([left1 +nGhostSkip*width1 ,left2+nGhostSkip*width2],axis=0)
    OverlapRight = na.amin([right1-nGhostSkip*width1,right2-nGhostSkip*width2],axis=0)


    if debug > 1:
        print "Overlap left (%6f,%6f,%6f) right (%6f,%6f,%6f)"%\
              tuple(OverlapLeft.tolist()+OverlapRight.tolist())        
    start1 = nint( (OverlapLeft - left1)/width1)
    end1 =   nint( (OverlapRight- left1)/width1) #the -1 is taken care of in the slice object
    start2 = nint( (OverlapLeft - left2)/width2)
    end2 =   nint( (OverlapRight- left2)/width2) 


    if debug > 1:
        print "Out of the box"
        print "start1 (%4d,%4d,%4d) end1 (%4d,%4d,%4d)"%\
              tuple(start1.tolist()+end1.tolist())
        print "start2 (%4d,%4d,%4d) end2 (%4d,%4d,%4d)"%\
              tuple(start2.tolist()+end2.tolist())
    end1 = na.amax( [end1, start1 + 1.0] ,axis=0)
    end2 = na.amax( [end2, start2 + 1.0] ,axis=0)
    if debug > 1:
        print "After the bump"
        print "start1 (%4d,%4d,%4d) end1 (%4d,%4d,%4d)"%\
              tuple(start1.tolist()+end1.tolist())
        print "start2 (%4d,%4d,%4d) end2 (%4d,%4d,%4d)"%\
              tuple(start2.tolist()+end2.tolist())

    slice1 = [slice(start1[i],end1[i]) for i in range(start1.shape[0])]
    slice2 = [slice(start2[i],end2[i]) for i in range(start2.shape[0])]
    #print "CLOWN", slice1
    #print "CLOWN", slice2

    return (slice1,slice2)




