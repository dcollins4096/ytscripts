

import numpy as na
nar = na.array

from yt.mods import *
import copy
import FindOverlap
import types
reload(FindOverlap)
from davetools import *
import re
import h5py
from dsdiff_helpers import *

    
class udiff():

    def togglev(self):
        """toggle velocity"""
        toggle(self.p['fields'],['x-velocity','y-velocity','z-velocity'])
    def togglebc(self):
        """toggle centered magnetic fields"""
        toggle(self.p['fields'],['MagneticField_C_1','MagneticField_C_2','MagneticField_C_3'])
    def togglebf(self):
        """toggle centered magnetic fields"""
        toggle(self.p['fields'],['MagneticField_F_1','MagneticField_F_2','MagneticField_F_3'])
    def togglee(self):
        """toggle electric fields"""
        toggle(self.p['fields'],['ElectricField_1','ElectricField_2','ElectricField_3'])
    def togglea(self):
        """toggle accelerations"""
        toggle(self.p['fields'],['Accel0','Accel1','Accel2'])

    def __init__(self,dir1,dir2=None,sticky=True,**kwargs):
        """Set up comparisons between *dir1* and *dir2*.  If *uber2* is
        ommited, *uber1* is used for both.
        *sticky* causes arguments to the __call__ method to replace existing
        methods.  If false, the move is only temporary"""
        self.dir1 = dir1
        if dir2:
            self.dir2=dir2
        else:
            self.dir2=dir1
        self.output_format ="%s_n%04d_g%04d_%s_%s%02d_%s.png"

        self.p={}
        self.p['diff_extents']=None
        self.p['sticky']=sticky
        self.p['frames']=kwargs['frames']
        self.p['grids']=[1]
        self.p['fields']=kwargs['fields']
        self.p['GhostInImage']=True
        self.p['nGhost']=0
        self.p['shift']=nar([0.0, 0.0, 0.0])
        self.p['output']=None
        self.p['raster']=None
        self.p['interogate_range']=[0]
        self.p['interogate_at_max']=True
        self.p['interogate_at_half']=True
        self.p['normalize']=True
        self.p['output_prefix']='dbmf'
        self.p['grid_direct']=True

        self.update(kwargs)
    def __getitem__(self,item):
        return self.p[item]
    def __setitem__(self,item,value):
        self.p[item] = value
    def update(self,kwargs):
        """Compare kwargs to the current dict.
        If self.p (the parameter set) already has the key, update.
        Otherwise, issue warning and proceed as normal."""
        pkeys = self.p.keys()
        for key in kwargs.keys():
            if key in pkeys:
                self.p[key] = kwargs[key]
            else:
                print "Unknown keyword: (%s) ignoring."%key
#       for f in ['frames','grids','fields']:
#           self.p[f] = ensure_list(self.p[f])
    def __call__(self,**kwargs):
        """
        kwargs can be:
        sticky      : if True, all kwargs are set permanently
        frames, grids: lists of frames and grids for each uber.
                      if an element in the list is a tuple, then the first
                      element of that tuple is used for uber1, the second for uber2.
                      e.g. grids=[(1,2),3]
                      will compare uber1/grid 1 to uber2,grid 2
                      then will compare uber1/grid3 to uber2/grid3
        fields      : list of fields to compare.  See the toggle methods for easy swapping
        nGhost=0    : Number of ghost zones in the grids
        GhostInImage: If False, *nGhost* zones are stripped from the image
        shift       : amount to shif along each axis. default = [0,0,0]
        output      : =None,'x','y','z'.  Slices along that axis using *interogate_range*
                      and plots an image of both grids an the difference.
        raster      : =None,x,y,z.  Slices along that axis using *interogate_range*
                      and querries differences along each slice
        interogate_range : the range for output and raster slices.
        normalize   : use relative differences (a-b)/(0.5 (a+b))
        output_prefix : prefix for differences
        output_format : format for the output filenames.  Must accept
                        output_prefix, field,output,xyz,n1,grid1,uber1.outname
        """

        total = 0.0
        maxdiff = None
        max_raster=0.0

        skipping = ""

        for p in self.p.keys():
            """Set the value in the local namespace,
            using the value in kwargs first, self.p second"""
            string = "%s = kwargs.get('%s',self.p['%s'])"%(p,p,p)
            exec(string)
#       for field in ['frames','grids','fields']:
#           exec("%s = ensure_list(%s)"%(field,field))

        if self['sticky']:
            self.update(kwargs)
        for n in frames:
            try:
                n1 = n[0]
                n2 = n[1]
            except:
                n1=n
                n2=n
            ds_name_1 = get_ds_name(self.dir1, n1)
            ds_name_2 = get_ds_name(self.dir2, n2)
            if not self.p['grid_direct']:
                ds1 = yt.load(ds_name_1)
                ds2 = yt.load(ds_name_2)
            for g in grids:
                try:
                    grid1 = g[0]
                    grid2 = g[1]
                except:
                    grid1=g
                    grid2=g
                for field in fields:
                    print "=================",field,"================="
                    if isinstance(field, types.ListType):
                        field1 = field[0]
                        field2 = field[1]
                    else:
                        field1=field
                        field2=field
                    try:
                        #Slice1, Slice2 = FindOverlap.FindOverlap(self.uber1,self.uber2,n1,n2,grid1,grid2,field1,shift,nGhost,self.p['grid_direct'])
                        Slice1 = slice(None)
                        Slice2 = slice(None)
                    except FindOverlap.OverlapException as e:
                        print e
                        continue
                    except Exception as e:
                        raise
                    if self.p['grid_direct']:
                        g1_full=  read_grid(self.dir1,n1,grid1,field1)
                        g2_full = read_grid(self.dir2,n2,grid2,field2)
                    else:
                        #verbose_save = self.uber1.verbose, self.uber2.verbose
                        #self.uber1.verbose, self.uber2.verbose=False,False
                        #self.uber1.fill(n1,get_region=False)
                        #self.uber2.fill(n2,get_region=False)
                        g1_full = self.uber1.h.grids[grid1-1][field1]
                        g2_full = self.uber2.h.grids[grid2-1][field2]
                        #self.uber1.verbose, self.uber2.verbose = verbose_save 
                    g1 = g1_full[Slice1]
                    g2 = g2_full[Slice2]

                    self.g1_full=g1_full; self.g2_full=g2_full
                    self.g1=g1; self.g2=g2
                    self.Slice1=Slice1; self.Slice2=Slice2

                    if g1.shape != g2.shape:
                        print "Shape difference: g1", g1.shape, "g2", g2.shape
                        print "Shape difference: g1", g1.shape, "g2", g2.shape
                        print "Shape difference: g1", g1.shape, "g2", g2.shape
                        print "Shape difference: g1", g1.shape, "g2", g2.shape
                        print "Shape difference: g1", g1.shape, "g2", g2.shape
                        continue

                    self.diff = (g1-g2)
                    diff = self.diff
                    if self.p['normalize']:
                        denom = (g1+g2)
                        ok = denom != 0
                        denom = denom[ ok ]
                        diff[ ok ] /= (0.5*(denom))
                    thismax = na.abs(diff).max()
                    total += thismax

                    if maxdiff == None or maxdiff[0] < thismax:
                        maxdiff = (thismax,n1,n2,field1,field2,grid1,grid2)

                    stat(g1,"g1: n1 = %d g1 = %d f = %s %s"%(n1,grid1,field1,'set1'))
                    stat(g2,"g2: n2 = %d g2 = %d f = %s %s"%(n2,grid2,field2,'set2'))
                    stat(diff,"n = %s g = %s f = %s, %s"%(str(n),str(g), field1, field2))

                    if raster:
                        max_this=0
                        for x in interogate_range:
                            subset = [slice(None)]*3
                            index = {'x':0,'y':1,'z':2}[raster]
                            subset[index]=x
                            d2=diff[subset]
                            stat(d2,"   %s %d"%(raster,x))
                            max_raster = max(max_raster,na.abs(d2).max())
                            max_this = max(max_this,na.abs(d2).max())
                        print "max_raster_this",max_this

                    if output:
                        if self.p['diff_extents'] is not None:
                            diff = na.clip(diff, self.p['diff_extents'][0], self.p['diff_extents'][1])
                        if not GhostInImage:
                            slT = slice(nGhost,-nGhost)
                        else:
                            slT = slice(None)
                        subset = [slT]*3
                        index = {'x':0,'y':1,'z':2}[output]
                        if self.p['interogate_at_max']:
                            all_max = na.array(na.where(na.abs(self.diff)==na.abs(self.diff).max()))
                            first_max = [stripe[0] for stripe in all_max]
                            interogate_range=[first_max[index]]
                            print "Peak difference at ", interogate_range
                        if self.p['interogate_at_half']:
                            interogate_range=[self.diff.shape[index]/2]
                            print interogate_range
                        for x in interogate_range: #range(g1.shape[0]):
                            subset[index]=min([x,g1[subset].shape[index]-1,g2[subset].shape[index]-1])
                            plave(g1[subset],self.output_format%(output_prefix,\
                                                           n1,grid1,field1,output,x,'set1'))
                            plave(g2[subset],self.output_format%(output_prefix,
                                                           n2,grid2,field2,output,x,'set2'))
                            diff_out = field1
                            if field1 != field2:
                                diff_out += field2
                            plave((diff)[subset],self.output_format%(output_prefix,\
                                                               n1,grid1,diff_out,output,x,'diff'))

                    print "total =",total, "max", maxdiff
                    if raster:
                        print "max_raster",max_raster
