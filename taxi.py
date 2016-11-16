"""The taxi object, for keeping track of yt things."""
""" used to be called uber, but that got taken. """
#checkrel
import warnings
not_ported = False
#with warnings.catch_warnings():
    #warnings.simplefilter("ignore")
    #import yt.raven as raven
    #import yt.lagos as lagos
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
import numpy as np
import pdb
import os
import h5py
import types, time,weakref,davetools
import dave_callbacks
from davetools import dsave, no_trailing_comments, ensure_list
import time, fPickle, glob
if not_ported:
    import clump_stuff, clump_subset
    from clump_list_sort import *
    #reload(clump_stuff)
import matplotlib.colorbar as cb
import re
import copy

def lim_down(value):
    return 10**(np.floor(np.log10(value)))
def lim_up(value):
    return 10**(np.ceil(np.log10(value)))

class PortException(Exception):
    def __init__(self, value = ""):
        self.value = "Taxi function not yet ported" + value
    def __str__(self):
        return repr(self.value)
class OperationException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

#GET def callbacks():
    
#GET def operations():

def writefunction(thing):
    """
    uber.read(file) populates itself from file by executing each line in the file.
    writefunction(thing) returns a string that, when executed, repopulates uber.
    Assumes all numpy nd arrays are 1d.  This can be fixed with recursion, but I'm lazy."""

    output = ""
    if isinstance(thing, types.ListType):
        output += "["
        for L in thing:
            output += writefunction(L)
            output += ", "

        #There are probably too many commas.
        if len(thing) > 0:
            output = output[0:-2]
        output += "]"
    elif isinstance(thing, types.StringType):
        output += "'"+thing+"'"
    elif isinstance(thing, np.ndarray):
        tmp = "np.array(["
        for L in thing:
            tmp += str(L)
            tmp += ", "
        output += tmp[0:-2] + "])"
    else:
        output += str(thing)
    return output

class taxi:
    """Taxi is a container for yt options, because the author is both forgetful and lazy.
    Primarily use for repeatability and reduced startup for plotting simple things that I do
    all the time, with default values suitably chosen.  See the source code for the options
    and defaults.  The won't be repeated here because there are several, and that's asking
    for a documentation inconsistency."""
    
    def __init__(self, filename=None,dir=None, name=None,**kwargs):
        """
        Either reads in taxifile *fileame* or just sets some defaults."""

        #List of members that will NOT be written.  This is mostly for the lagos objects.
        #The goal is to be able to re-create the lagos objects from the info in uber,
        #and not have to pickel the lagos objects.
        self.ExcludeFromWrite = ['ExcludeFromWrite']

        #To make the file easier to use, some of the uber members are written first
        self.WriteMeFirst = ['name','directory','outname','operation','frames','fields','axis',
                            'name_syntax', 'setname','DirectoryPrefix','GlobalParameters',
                             'ProfileDir','ProfileName']
        self.ExcludeFromWrite.append('WriteMeFirst')

        #Some members may need special treatment.
        self.WriteSpecial = {}
        self.ExcludeFromWrite.append('WriteSpecial')
        #externals
        self.OutputLog = "OutputLog"
        self.plot_origin = 'domain'  #see yt docs.
        self.axis = 0                #Axis for the plot
        self.frames = []             #Dump numbers to be plotted.
        self.fields = 'Density'      #Field (hopefully to be list) for the plot
        self.weight = 'dx'           #Weight field.  Projections?    
        self.name   = 'sim'          #identifier for this uber instance
        self.outname = 'Image'       #Prefix for output
        self.directory = '.'         #Directory where the simulation was run. (Where the DD* dirs are)
        self.DirectoryPrefix = "DD"  #Prefix on said directories
        self.ProfileDir= "./ProfileFiles/"
        self.ProfileName = None
        self.setname = "cycle"        #prefix on data file
        self.subdir = True           #Are there DD0000 directories, or is it just data0000.grid*?
        self.operation = 'Full'      #The operation that will be done.  uber.operations() for a list.
                                     #    or physically motivated ones.  Options are 'Code' or 'Physics'
        self.format = 'png'
        self.plot_args = {}
        #The actual yt objects.  Not saved on output.
        self.ds = None               #The lagos parameter file.  Gets re-made at each plot.
        self.ExcludeFromWrite.append('ds')

        self.pc = None               #The plot collection.
        self.ExcludeFromWrite.append('pc')
        self.index = None                #A weak proxy to ds.h
        self.ExcludeFromWrite.append('index')
        self.h = None                #A weak proxy to ds.h
        self.ExcludeFromWrite.append('h')
        self.reg = None              #Weak proxy to the most recent region
        self.ExcludeFromWrite.append('reg')
        self.proj = None             #Weak proxy to the most recent projection
        self.ExcludeFromWrite.append('proj')

        self.extra_save = False       #Forces an extra plot save to get colormaps correct.
        self.ExcludeFromWrite.append('extra_save')        

        self.ExcludeFromWrite.append('region')
        self.field_parameters = {}
        self.left=np.array([0.0,0.0,0.0])      #The left of the region extraction.
        self.right=np.array([1.0,1.0,1.0])     #The right of the extracted region.
        self.radius=None
        self.restrict = False                   #Restric plot to img_left, img_right
        self.img_left = None                   #left edge of image
        self.img_right = None                  #right edge of image.
        self.img_width = None                  #width of image.
        self.width = None
        self.region_type = 'sphere'
        self.center= 0.5*(self.left+self.right)#Center, used for all things that need a center
        self.periodic = False              #For projection shift
        self.cmap = {'default':None}       #color map, field driven. Everyone uses default.
        self.Colorbar = None               #Monotone: monotonically increaseing, forced to min,max
        self.use_colorbar = True           #use color bar for plots.  Bad naming.
        self.zlim = {}
        self.slice_zlim = {}
        self.proj_zlim = {}
        self.callbacks = []                #Which callbacks to use.  uber.callbacks() for options
        self.callback_args={}
        #self.ExcludeFromWrite.append('callbacks')
        self.WriteSpecial['callbacks'] = self.write_callbacks

        self.MaxLevel=-1                   #max analysis level, currently does nothing
        
        #for strings on the image
        self.UnitMultiplier = 1.0          #Info for the PrintTime callback.  Kind of a hack.
        self.UnitName = 'FemtoParsec'      #Info for the PrintTime callback.  Kind of a hack.
        
        #for clump reading
        self.clump_prefix = None
        self.clumps = None
        self.ExcludeFromWrite.append('clumps')
        
        #periodic?
        self.Periodic = False              #For periodically shifting boxes when plotting.
        
        #internals
        self.old_outname = None
        self.frame_template = None
        self.basename_template = "data"
        self.name_syntax = 'outputlog'
        self.extrema = {}
        self.frame_dict = None
        self.ExcludeFromWrite.append('frame_dict')
        self.basename = None
        self.verbose = True
        self.mark_time = True
        #for storing global parameters
        #Updates each parameter file on read.
        self.GlobalParameters = {}

        #Covering grid.
        self.cg = None
        self.ExcludeFromWrite.append('cg')
        self.hold_values={}
 


        if filename != None:
            #If no directories are mentioned, assume the actual file is ./uber/filename
            try:
                filename.index('/')
                actual_filename = filename
            except:
                actual_filename = "uber/%s"%filename
            self.load(actual_filename)
        if dir != None:
            self.directory=dir
        if name != None:
            self.name = name
            self.outname = name
        self.__dict__.update(**kwargs)

    def load(self,filename='DefaultFile'):
        """Simply executes each line in the file *filename* in the namespace of this uber."""

        file = open(filename,"r")
        for line in file:
            exec(line)
        file.close()

    def save(self,filename='DefaultFile'):
        """Writes out the taxi file.  Writes the entire contents of taxi.__dict__,
        making a hopefully appropriate assumption about the type."""
        file = open(filename,"w")
        for i in self.WriteMeFirst:
            file.write("self."+i + " = " + writefunction(self.__dict__[i]) + "\n")            

        for i in self.__dict__.keys():
            if i in self.ExcludeFromWrite:
                continue
            if i in self.WriteMeFirst:
                continue
            if i in self.WriteSpecial.keys():
                continue
            file.write("self."+i + " = " + writefunction(self.__dict__[i]) + "\n")
        for i in self.WriteSpecial.keys():
            file.write( self.WriteSpecial[i]() )
        file.close()

    def write_callbacks(self):
        """checkrel"""
        output = []
        for i in self.callbacks:
            if isinstance(i,types.StringType):
                output.append(i)
        return "self.callbacks = " + writefunction(output)


    def __str__(self):
        out = ""
        for i in self.WriteMeFirst:
            out += "self."+i + " = " + writefunction(self.__dict__[i]) + "\n"
        for i in self.__dict__.keys():
            if i not in self.ExcludeFromWrite and i not in self.WriteMeFirst:
                out += "self."+i + " = " + writefunction(self.__dict__[i]) + "\n"
        return out

    # checkrel def frame_grabber(self,tinitial=None,tfinal=None,dt=None,preset=None):

    #
    #
    # filenames are a mild hassle.
    # Either they can be preset, e.g DD????/data???? where DD and data are variables.
    # OR can be taken sequentially from the OutputLog.
    # 

    def set_filename(self,frame=None):
        if self.name_syntax == 'preset':
            self.set_filename_preset_syntax(frame)
        elif self.name_syntax == 'outputlog':
            self.set_filename_from_outputlog(frame)
        else:
            print "set oober.name_syntax = 'preset'||'outputlog'"


    def return_filename(self,frame=None):
        verb_save=self.verbose
        self.verbose=False
        self.set_filename(frame)
        self.berbose=verb_save
        return self.basename

    def set_filename_preset_syntax(self,frame=None):
        if frame == None:
            frame = self.frames[-1]
        #take 1: no subdirectory (very old style)
        self.basename_template = self.directory+'/'+self.setname+'%04i'
        self.basename = self.basename_template %(frame)
        if glob.glob(self.basename+".hierarchy") == []:
            self.basename_template = self.directory + "/" + self.DirectoryPrefix + '%04i/'\
                                     + self.setname + "%04i"
            self.basename = self.basename_template %(frame,frame)
        if self.verbose == True:
            print "==== %s ===="%(self.basename)

    def print_frames(self):
        self.get_frames()
        if self.name_syntax == 'preset':
            print "oober.name_syntax == 'preset'.  Nothing to list"
        elif self.name_syntax == 'outputlog':
            self.get_frames()
            nframe = self.frame_dict.keys()
            nframe.sort()
            output_format = "%s%d%s"%("%(dsname)",self.frame_dict_string_length ,"s %(cycle)d %(time)0.6e")

            for n in nframe:
                print "%4d"%n, output_format%(self.frame_dict[n])
        else:
            print "set oober.name_syntax = 'preset'||'outputlog'"

    def returnsetnumber(self,frame=None):
        return frame
        if self.name_syntax == 'preset':
            return frame
        else:
            self.get_frames()
            if frame==None:
                frame = self.frames[0]
            return self.frame_dict[frame]['SetNumber']

    def get_frames(self):
        """parses OutputLog to populate self.frame_dict"""
        DirSetNumber = re.compile(r'([^\d]*)(\d\d\d\d)/([^\d]*)(\d\d\d\d)')
        SetNumber = re.compile(r'([^\d]*)(\d\d\d\d)')
        SetNumber = re.compile(r'(.*)(\d\d\d\d)')
        debug = 0
        if not hasattr(self,'frame_dict') or self.frame_dict is None:
            self.frame_dict = {}
            output_log = open("%s/%s"%(self.directory , self.OutputLog),'r')
            if debug > 0:
                print output_log
            lines = output_log.readlines()
            self.frame_dict_string_length = 0
            for nframe,line in enumerate(lines):
                if debug>0:
                    print nframe, line[:-1]
                linesplit = davetools.no_whites(line.split(" "))
                self.frame_dict[ nframe ] = {}
                self.frame_dict[ nframe ]['dsname'] = linesplit[2]
                self.frame_dict_string_length = max( self.frame_dict_string_length , len(linesplit[2]) )

                try:
                    self.frame_dict[ nframe ]['cycle']  = int(linesplit[3])
                except:
                    self.frame_dict[ nframe ]['cycle']  = -1
                match = DirSetNumber.match(linesplit[2])
                match2 = SetNumber.match(linesplit[2])
                if match is not None:
                    self.frame_dict[nframe]['DirPrefix']=match.group(1)
                    self.frame_dict[nframe]['SetName'] = match.group(3)
                    self.frame_dict[nframe]['SetNumber'] = int(match.group(2))
                elif match2 is not None:
                    self.frame_dict[nframe]['DirPrefix']=None
                    self.frame_dict[nframe]['SetName'] = match2.group(1)
                    self.frame_dict[nframe]['SetNumber'] = int(match2.group(2))
                else:
                    print "error in set parsing"
                    print linesplit[2]

                    self.frame_dict[nframe]['DirPrefix']=None
                    self.frame_dict[nframe]['SetName']  =None
                    self.frame_dict[nframe]['SetNumber']=-1

                try:
                    self.frame_dict[ nframe ]['time']   = float(linesplit[4])
                except:
                    self.frame_dict[ nframe ]['time']   = -1
                #pdb.set_trace()
    def set_filename_from_outputlog(self,frame=None):
        self.get_frames()
        if frame == None:
            frame = self.frames[-1]
        self.basename = "%s/%s"%(self.directory,self.frame_dict[frame]['dsname'])
        if self.verbose == True:
            print "==== %s ===="%(self.basename)

    def fill(self, frame = None):
        """populate parameter file, plot collection, hierarchy, region if desired.
        Possibly could be renamed 'load' """
        self.set_filename(frame)
        self.ds = yt.load(self.basename)
        na_errors= np.seterr(all='ignore')
        np.seterr(**na_errors)
    def get_region(self, frame=None):
        if frame is not None or self.ds is None:
            self.fill(frame)
        if self.region_type.lower()=='sphere':
            reg = self.ds.sphere(self.center, self.radius)
        if self.region_type.lower() in ['rectangle','region']:
            self.center = 0.5*(self.left + self.right)
            reg = self.region(self.center, self.left,self.right,field)
        if self.region_type.lower() in ['all','all_data']:
            reg = self.ds.all_data()
        #self.reg  = weakref.proxy(reg)
        return reg




    def plot(self,local_frames=None):
        print self.zlim
        """Makes the uber plot.  Options for plot can be found in the uber description,
        plot operations can be found by typing uber.operations()"""
        start_time = time.time()    
        
        frame_template = self.outname + "_%04i"
        FirstOne = True #For monotonic plots.

        if local_frames==None:
            local_frames=self.frames

        for frame in ensure_list(local_frames):

            self.fill(frame=frame)

            print "==== %s ===="%(self.basename)
            print "== %i time elapsed %f =="%(frame, time.time()-start_time)

            #Find the density peak for all field,axis plots for this snapshot.
            if self.operation == 'DensityPeakSlice':
                print "getting peak"
                value,center = self.ds.find_max('Density')
                self.center = np.array(center)

            for field in ensure_list(self.fields):

                #Find the center, if the center is to come from the min/max.
                #If a multi plot, use the first plot given.  
                
                if self.operation == 'PeakSlice':
                    value,center = self.ds.find_max( ensure_list(field)[0] )
                    self.center = np.array(center)
                elif self.operation == 'MinSlice':
                    value,self.center = self.ds.find_min( ensure_list(field)[0] )
                    self.center = np.array(center)

                #start axing
                for axis in ensure_list(self.axis):
                    #if multiple plots are used, do_plot will return a figure object.
                    plot_or_fig = self.do_plots(field,axis,FirstOne)
                    if hasattr(plot_or_fig,'savefig'):
                        
                        filename = frame_template%(self.returnsetnumber(frame))
                        filename += "_"+self.operation
                        filename += "_"+['x','y','z'][axis]

                        for i in field:
                            filename += "_"+i
                        if self.format:
                            filename += "."+self.format
                        else:
                            filename += ".png"
                        print filename
                        #fig.savefig(filename)
                        #dsave(fig,filename,field_name=field,ds_list=[self.basename],script_name='uber')
                        plot_or_fig.savefig(filename)
                        print filename
                    else:
                        print plot_or_fig.save(frame_template%frame)
    def do_plots(self,sub_field,axis,FirstOne):

        """Adds the plot to the plot collection.  Loops over plots for multi-plots.
        FirstOne tracks the first plot, for colorbar purposes"""
        sub_field_list = ensure_list(sub_field)
        orient = 'horizontal' #hard coded for devel!  
        use_colorbar = self.use_colorbar #local copy made for multi plots
        extra_args = self.plot_args

        #set up the multi plot if this field is a list.
        if len(sub_field_list) > 1:
            config = guess_multi_configuration(len(sub_field_list))
            #fig, axes, colorbars = raven.get_multi_plot( config[0], config[1], colorbar=orient, bw = 4)
            fig, axes, colorbars = get_multi_plot( 2, 1, colorbar=orient, bw = 4)
            use_colorbar = False 
            extra_args = {'figure':fig,'axes':0}
        for num, field in enumerate(sub_field_list):
            if extra_args.has_key('axes'):
                #extra_args['axes'] = axes[wp(num)[0]][wp(num)[1]]
                extra_args['axes'] = wp(axes,num) #axes[wp(num)[0]][wp(num)[1]]

            if self.operation == 'multi_test':
                the_plot = yt.SlicePlot(self.ds,axis,field,center=self.center,origin=self.plot_origin,**extra_args)
                #not tested.
                #the_plot = self.pc.add_slice(field,axis, use_colorbar=use_colorbar,
                #                             **extra_args)
            elif self.operation == 'Full':
                if not extra_args.has_key('weight_field'):
                    extra_args['weight_field'] = None
                the_plot = yt.ProjectionPlot(self.ds,axis,field, center=self.center, **extra_args)
                #self.pc.add_projection(field,axis, use_colorbar=use_colorbar,
                #                            **extra_args)
                self.zlim = self.proj_zlim

            elif self.operation == 'RegionProjection_kill_this':
                """Might be broken"""
                #setup a projection
                self.center = 0.5*(self.left + self.right)
                reg = self.region(self.center, self.left,self.right,field)
                self.reg = weakref.proxy(reg)
                self.proj = self.ds.proj(axis,field,center=self.center,source=self.reg,periodic=self.periodic)
                the_plot = self.proj.to_pw()
                self.zlim = self.proj_zlim
            elif self.operation == 'RegionProjection':
                """Might be broken"""
                #setup a projection
                reg = self.get_region()
                self.proj = self.ds.proj(field,axis,center=self.center,data_source=reg) #,periodic=self.periodic)
                the_plot = self.proj.to_pw()
                self.zlim = self.proj_zlim
            elif self.operation in ['MinSlice','PeakSlice','DensityPeakSlice','CenterSlice']:
                """peak taken outside the axis loop"""
                the_plot = yt.SlicePlot(self.ds,axis,field,center=self.center,origin=self.plot_origin,**extra_args)
                self.zlim = self.slice_zlim
            else:
                raise OperationException(self.operation)

            #this is to clean up old implementations
            if not self.cmap:
                self.cmap = {'default':None}
            if not self.cmap.has_key(field):
                self.cmap[field] = self.cmap['default']

            the_plot.set_cmap( field, self.cmap[field] )
            #the_plot.label_kws['size'] = 'x-large'
            if self.Colorbar:
                if self.Colorbar in ['Monotone', 'monotonic']:
                    if FirstOne:
                        self.zlim[field] = [the_plot.data_source[field].min(),the_plot.data_source[field].max()]
                    else:
                        self.zlim[field][0] = min([self.zlim[field][0],the_plot.data_source[field].min()])
                        self.zlim[field][1] = max([self.zlim[field][1],the_plot.data_source[field].max()])
                    if self.verbose:
                        print "set lim", self.zlim[field]
                elif self.Colorbar == 'Fixed':
                    if not self.zlim.has_key(field):
                        self.zlim[field] = [the_plot.data_source[field].min(),the_plot.data_source[field].max()]
                    if self.verbose:
                        print "set lim", self.zlim[field]
                the_plot.set_zlim(field, self.zlim[field][0], self.zlim[field][1])

                if self.operation in ['Full','RegionProjection']:
                    self.proj_zlim = self.zlim
                elif self.operation in ['MinSlice','PeakSlice','DensityPeakSlice','CenterSlice']:
                    self.slice_zlim = self.zlim

        self.the_plot = the_plot #this does actually need to be a list.
        FirstOne = False
        try:
            self.add_callbacks(the_plot, axis=axis)
        except:
            pass
        if self.restrict:
            if self.width:
                the_plot.set_width(self.width)                    

        if len(sub_field_list) > 1:
            raise PortException("multiplots")
            for p,cax in zip(self.pc.plots, colorbars):
                # Now we make a colorbar, using the 'image' attribute of the plot.
                # 'image' is usually not accessed; we're making a special exception here,
                # though.  'image' will tell the colorbar what the limits of the data are.
                cbar = cb.Colorbar(cax, p.image, orientation=orient)
                # Now, we have to do a tiny bit of magic -- we tell the plot what its
                # colorbar is, and then we tell the plot to set the label of that colorbar.
                p.colorbar = cbar
                p._autoset_label()                
            return fig
            
        else:
            return the_plot

    def product_dir(self,frame):
        """Destination of intermediate products.  Typically of the form 
        simdir/DD0001.products
        """
        if self.name_syntax == 'preset':
            return self.pdir_old(frame)
        elif self.name_syntax == 'outputlog':
            self.get_frames()
            return "%s/%s%04d.products"%(self.directory,
                    self.frame_dict[frame]['DirPrefix'],
                    self.frame_dict[frame]['SetNumber'])

    def pdir_old(self,frame):
        """Finds .products directores in self.directory.  
        Makes the questionable assumption that directories are uniquely numbered."""
        print "old style pdir needs to be checked."
        raise PortException("old style product directory")
        if hasattr(self,'pdirdict'):
            return self.pdirdict[frame]
        else:
            produs = glob.glob(self.directory+"/*.products")
            first_to_check = "%s/%s%04d.products"%(self.directory,
                                                   self.DirectoryPrefix,
                                                   frame)

            if first_to_check in produs:
                return first_to_check

            self.pdirdict = {}
            TheNumbers = re.compile(r'.*(\d\d\d\d).products')
            for L in produs:
                pdir = L.split("/")[-1]
                match = TheNumbers.match(pdir)
                if match is not None:
                    self.pdirdict[ int(match.groups(1)[0])] = L
            return self.pdirdict[frame]

    def make_cg(self,frame=None,field=None, num_ghost_zones=0):
        """Makes a covering grid for the root grid."""
        if field == None:
            field = self.fields[0]
        self.fill(frame)
        left=[0.0]*3
        resolution = self.ds['TopGridDimensions'] + 2*num_ghost_zones
        self.cg=self.ds.covering_grid(0,left,resolution,num_ghost_zones=num_ghost_zones)#,fields=fields)

    def extract_root(self,frame=None,field=None,make_cg=True,num_ghost_zones=0,debug=0,dtype='float32'):
        """Extracts a grid the size of the root grid.  
        Cached in self.directory+".products"+"cube_<field>.single" exists, reads that.
        *num_ghost_zones* is for fields that require derivatives, etc.
        Currently defaults to first frame, all fields"""
        if frame == None:
            frame = self.frames[0]
        if field == None:
            field = self.fields[0]
        directory = self.product_dir(frame)
        filename_no_ghost = "%s/cube_%s.%s"%(directory,field,dtype)
        if debug > 0:
            print filename_no_ghost
        remove_gz = False
        if num_ghost_zones >= 0:
            filename = "%s.%d"%(filename_no_ghost,num_ghost_zones)
        else:
            """If num_ghost_zones < 0, try to read with ngz=0.  If not there, try larger values.
            If extant, use that and remove ghost zones."""
            filename = "%s.%d"%(filename_no_ghost,0)
            if glob.glob(filename) == []:
                with_ghost_list = glob.glob(filename_no_ghost+"*")
                if with_ghost_list != []:
                    ghost_zones_available = [ int(a.split(".")[-1]) for a in with_ghost_list]
                    ngz = min(ghost_zones_available)
                    filename = "%s.%d"%(filename_no_ghost, ngz)
                    if debug>0:
                        print "Found",filename, "with ngz",ngz
                    remove_gz = True

        if glob.glob(filename):
            if debug > 0:
                print "open set from disk"
            file = h5py.File(filename,'r')
            set = file[field][:]
            file.close()
            if remove_gz:
                if debug>0:
                    print "stripping"
                set = set[ngz:-ngz,ngz:-ngz,ngz:-ngz]
        else:
            if glob.glob(directory) == []:
                print "making directory",directory
                os.mkdir(directory)
            if debug > 0:
                print "Create set"
            if make_cg or self.cg == None:
                if debug > 0:
                    print "make cg"
                self.make_cg(frame,field,num_ghost_zones)
            set = self.cg[field]
            file = h5py.File(filename,'w')
            file.create_dataset(field,set.shape, data=set, dtype=dtype)
            file.close()
        return set

    def fft(self,frame=None,field=None,data=None,make_cg=True,num_ghost_zones=0,dtype='float32',debug=0,fft_func=np.fft.fftn):
        """Make an fft of a field.
        As with extract_root, look for a cached version on disk first.
        Defaults to the ghost zone guess"""
        if dtype == 'float32':
            fft_dtype = 'complex64'
        elif dtype == 'float64':
            fft_dtype = 'complex128'
        elif dtype not in ['complex64','complex128']:
            print "Can't cast type ",dtype, "to a complex type."
            return None
        if frame == None:
            frame = self.frames[0]
        if field == None:
            field = self.fields[0]
        directory = self.product_dir(frame) #"%s/%s%04d.products/"%(self.directory,self.DirectoryPrefix,frame)
        filename = "%s/fft_%s.%s"%(directory,field,dtype)
        if debug > 0:
            print filename
        if glob.glob(filename):
            if debug > 0:
                print "open FFT from disk"
            file = h5py.File(filename,'r')
            fft = file[field][:]
            file.close()
        else:
            if glob.glob(directory) == []:
                print "making directory",directory
                os.mkdir(directory)
            if debug > 0:
                print "Create FFT"
            if data == None:
                this_set = self.extract_root(frame,field,make_cg,num_ghost_zones,dtype)
            else:
                this_set = data
            if debug > 0:
                print "Do fft."
            fft = fft_func(this_set)/this_set.size
            if debug > 0:
                print "save"
            fptr = h5py.File(filename,'w')

            fptr.create_dataset(field,fft.shape, data=fft, dtype=fft_dtype)
            fptr.close()
        return fft
    def find_extrema(self,fields=None,frames=None):
        if fields:
            local_fields = fields
        else:
            local_fields = self.fields
        if frames:
            local_frames = frames
        else:
            local_frames = self.frames
        for frame in local_frames:
            self.fill(frame)
            reg=self.get_region(frame)
            for field in local_fields:
                this_extrema = reg.quantities['Extrema'](field)
                if not self.extrema.has_key(field):
                    self.extrema[field] = [this_extrema[0].v, this_extrema[1].v]
                else:
                    self.extrema[field][0] = min([this_extrema[0].v, self.extrema[field][0]])
                    self.extrema[field][1] = max([this_extrema[1].v, self.extrema[field][1]])

    def phase(self,fields, callbacks=None, weight_field=None, phase_args={},save=True, n_bins=[64,64]):
        """Uber wrapper for phase objects.
        for all frames in self.frame, run a phase plot object on the region.
        Run all callbacks on the plot.
        *save* pickles the region in PhaseFiles"""
        if len(fields) != 3:
            print "wrong number of fields: needs 3.", fields
            return
        frame_template = self.outname + "_%04i"
        for frame in ensure_list(self.frames):
            reg = self.get_region(frame)
            local_extrema = None
            if self.Colorbar in ['fixed', 'monotonic']:
                print "start: ", self.extrema
                if self.extrema.has_key(fields[0]) and self.extrema.has_key(fields[1]):
                    local_extream = self.extrema
                else:
                    local_extrema = None
            print "mid: ", local_extrema
            phase = yt.create_profile(reg,bin_fields=[fields[0],fields[1]], 
                                      fields=[fields[2]],weight_field=weight_field,
                                      extrema=local_extrema, n_bins=n_bins)
            #self.phase = weakref.proxy(phase)
            pp = yt.PhasePlot.from_profile(phase)
            pp.set_xlabel(fields[0])
            pp.set_ylabel(fields[1])
            if 1:
                if self.Colorbar == 'monotonic':
                    if self.extrema.has_key(fields[0]):
                        xmin=min([self.extrema[fields[0]][0], phase.x_bins[0]])
                        xmax=max([self.extrema[fields[0]][1], phase.x_bins[-1]])
                        ymin=min([self.extrema[fields[1]][0], phase.y_bins[0]])
                        ymax=max([self.extrema[fields[1]][1], phase.y_bins[-1]])
                    else:
                        xmin=phase.x_bins[0]
                        xmax=phase.x_bins[-1]
                        ymin=phase.y_bins[0]
                        ymax=phase.y_bins[-1]
                    self.extrema[fields[0]] = [xmin,xmax]
                    self.extrema[fields[1]] = [ymin,ymax]
                    print "extrema, post", self.extrema
                if self.Colorbar in ['fixed','monotonic']:
                    pp.set_xlim( lim_down(self.extrema[fields[0]][0]), lim_up(self.extrema[fields[0]][1]))
                    pp.set_ylim( lim_down(self.extrema[fields[1]][0]), lim_up(self.extrema[fields[1]][1]))

            print pp.save("%s_%04d"%(self.outname,frame))
            old_code = """
            phase = self.pc.add_phase_object(self.region,fields, lazy_reader=True,weight=weight,**phase_args)
            if callbacks:
                for call in ensure_list(callbacks):
                    phase.add_callback(call)
                    
            if save and davetools.ImRoot():
                if glob.glob(self.ProfileDir) == []:
                    os.mkdir(self.ProfileDir)
                    print "made directory", self.ProfileDir
                filename = "%s/%s_%s_%s_%s_%s"%(self.ProfileDir,frame_template%(self.returnsetnumber(frame)),fields[0],fields[1],fields[2],weight)
                all_files = glob.glob(filename+"*")
                if 0:
                    filename += "_%d.pickle"%(len(all_files))
                    #fPickle.bdump(self.pc.plots[-1].data._data,filename)
                    fPickle.bdump(self.pc.plots[-1].data,filename)
                    print "saved pickle %s"%filename
                if 1:
                    filename += "_%d.h5"%(len(all_files))
                    print "Writing H5 file",filename
                    fptr = h5py.File(filename,"w")
                    the_dict = self.pc.plots[0].data.field_data
                    for fld in the_dict.keys():
                        this_field = the_dict[fld]
                        fptr.create_dataset(fld,this_field.shape, data=this_field)
                    fptr.create_dataset('x_bin_field',data=self.pc.plots[0].data.x_bin_field)
                    fptr.create_dataset('y_bin_field',data=self.pc.plots[0].data.y_bin_field)
                    fptr.create_dataset('ds_name',data=self.basename.split("/")[-1])
                    fptr.create_dataset('InitialTime',data=self.ds['InitialTime'])
                    fptr.create_dataset('InitialCycle',data=self.ds['InitialCycleNumber'])
                    fptr.close()


            if self.format == None:
                print self.pc.save(frame_template%(self.returnsetnumber(frame)))
            elif self.format == 'dave':
                dsave(self.pc,frame_template%(self.returnsetnumber(frame)),
                      field_name = "%s"*len(fields)%tuple(fields),
                      ds_list=[self.basename],
                      script_name='uber')
            else:
                print self.pc.save(frame_template%(self.returnsetnumber(frame)), format = self.format)
                #print frame_template%(frame)
                """
    def add_callbacks(self,the_plot,axis=None):
        for i,callback in enumerate(self.callbacks):
            if isinstance(callback,types.StringType):
                if callback == 'velocity':
                    the_plot.annotate_velocity()
                else:
                    raise PortError("Callback %s not supported"%callback)
class other_horsecrap():
######################### stuff not ported
    def set_region(self,frame=None):
        raise PortError('set_regin')
        self.region = self.region(0.5*(self.left+ self.right), self.left, self.right)
        for parameter in self.field_parameters.keys():
            self.region.set_field_parameter(parameter, self.field_parameters[parameter])

    def add_callbacks(self,axis=None):
        """Adds specified callbacks in the uberInstance.callbacks array to the last plot.
        Options are either callback shortcut string or instance of an appropriate callback
        object.  Acceptable shortcuts can be found in uber.callbacks()"""
                
        #
        # callbacks
        #
        if len(self.callbacks):
            raise PortError('callbacks')
        
        for i,callback in enumerate(self.callbacks):
            if isinstance(callback,types.StringType):
                if 'Subgrids' in self.callbacks:
                    for p in self.pc:
                        try:
                            args = self.callback_args['Subgrids']
                        except:
                            args = {'alpha':1.0}

                        p.annotate_grids()
                if 'Grids' in self.callbacks:
                    for p in self.pc:
                        p.annotate_grids()
                if 'Time' in self.callbacks:
                    for p in self.pc.plots:
                        p.modify['time'](unit_mult=self.UnitMultiplier, unit= self.UnitName)
                        #annotate_text([0.05,0.05],r'$t=0.5 t_{\rm{ff}}$',text_args={'color':(1.0,1.0,1.0,0.5),'fontsize':31})
                if 'Size' in self.callbacks:
                    print "NO SIZE CALLBACK FOR UBER"
                    continue
                    print "Nope. Define callback"
                    raise ZeroDivisionError
                    for p in self.pc.plots:
                        p.add_callback(dave_callback.SizeMarker())    
                if 'Zoom' in self.callbacks:
                    print "NO ZOOM CALLBACK FOR UBER"
                    continue
                    print "Nope. Define callback"
                    raise ZeroDivisionError
                    if self.width:
                        for p in self.pc.plots:
                            p.add_callback(dave_callback.BoxCallback(self.center,self.width))
                    else:
                        print 'width not defined.'
                if 'clumps' in self.callbacks:
                    print "NO CLUMPS CALLBACK FOR UBER"
                    continue
                    for p in self.pc.plots:
                        if not self.callback_args.has_key('clumps'):
                            return
                        clumps_to_use = self.callback_args['clumps']['clumps']
                        if self.callback_args['clumps'].has_key('annotate'):
                            annotate=True
                        else:
                            annotate=False
                        if self.callback_args['clumps'].has_key('sort'):
                            if self.callback_args['clumps']['sort']:
                                clumps_to_use = clump_list_sort(clumps_to_use)
                        #p.add_callback(davecallback.ClumpContourCallbackDave(clump_list_sort(clumps_to_use),annotate=annotate ))
                        p.add_callback(davecallback.ClumpContourCallbackDave(clumps_to_use,annotate=annotate ))
                if 'Particles' in self.callbacks:
                    for p in self.pc:
                        try:
                            args = self.callback_args['Particles']
                        except:
                            #(axis, width, p_size=1.0, col='k', stride=1.0
                            if not axis: axis = 0
                            args = {'width':1.0}
                        #args['axis'] = axis
                        p.annotate_particles( **args)
                if 'Velocity' in self.callbacks:
                    print "NO VELOCITy CALLBACK FOR UBER"
                    continue
                    for p in self.pc.plots:
                        try:
                            args = self.callback_args['Velocity']
                        except:
                            args={}
                        p.add_callback(VelocityCallback(**args))

            else:
                print "No Callback For You"

    def slicer(self,n_slices):
        """Make n_slices.
        ToDo: deal with different axes"""
        save={}
        save['outname']=self.outname
        prefix = self.outname
        save['center']=self.center
        save['operation'] = self.operation
        self.operation='CenterSlice'
        save['callbacks'] = self.callbacks
        base_callbacks=copy.copy(self.callbacks)
        axis = ensure_list(self.axis)[0] 
        axis_text = {0:'x',1:'y',2:'z'}[axis]
        try:
            for n in self.frames:
                for i,x in enumerate( np.arange(0.0,1.0, 1./n_slices) ):
                    self.center[axis] = x
                    self.outname = prefix + '_%04d'%i
                    text = '%s = %0.2f'%(axis_text,x)

                    #self.callbacks = base_callbacks + [davecallback.PrintText(text)]
                    self.plot(n)
        except:
            raise
        finally:
            self.__dict__.update(save)
    

    def quantities(self,quantity,*args,**kwargs):
        """Uber wrapper for derived quantities.
        For all frames in self.frame, run the given quantity on the region defined by (self.left, self.right).
        Returns a {frame:value} dictionary"""
        output = {}
        if kwargs.has_key('frames'):
            frames=kwargs.pop('frames')
        else:
            frames = self.frames

        for frame in frames:
            self.fill(frame)
            output[frame] = self.region.quantities[quantity](*args,**kwargs)
        return output


    def profiles(self,fields, callbacks=None, weight=None, profile_args={},
                 title=None,time=True,save=True,use_cg=False,lazy_reader=True,num_ghost_zones=0):
        """Uber wrapper for profile objects.
        for all frames in self.frame, run a profile object on the region.
        Run all callbacks on the plot.
        *use_cg* generates a covering grid and uses that, instead."""
        raise PortException('profiles')
        if len(fields) != 2:
            print "wrong number of fields: needs 2."
            return
        frame_template = self.outname + "_%04i"
        frame_template += "_%s"%weight
        for frame in ensure_list(self.frames):
            if use_cg:
                lazy_reader=False
                self.make_cg(frame=frame,num_ghost_zones=num_ghost_zones)
                obj=self.cg
                obj_string='_cg'
            else:
                self.fill(frame)
                obj=self.region
                obj_string=''

            if davetools.ImRoot():
                print 'ugly', fields, weight
            profile = self.pc.add_profile_object(obj,fields, lazy_reader=lazy_reader,weight=weight,**profile_args)
            if davetools.ImRoot():
                if self.mark_time:
                    mark_time("Made profile n=%d fields=%s,%s"%(self.returnsetnumber(frame),fields[0],fields[1]))
                if save:
                    if glob.glob(self.ProfileDir) == []:
                        os.mkdir(self.ProfileDir)
                        print "made directory", self.ProfileDir
                    filename = "%s/%s_%s_%s%s"%(self.ProfileDir,frame_template%(self.returnsetnumber(frame)),fields[0],fields[1],obj_string)
                    all_files = glob.glob(filename+"*")
                    filename += "_%d.h5"%(len(all_files))
                    file=h5py.File(filename,"w")
                    for i in range(2):
                        set = self.pc.plots[-1].data[fields[i]]
                        file.create_dataset(fields[i],set.shape,data=set)
                    file.create_dataset('ds_name',data=self.basename.split("/")[-1])
                    file.create_dataset('InitialTime',data=self.ds['InitialTime'])
                    file.create_dataset('InitialCycle',data=self.ds['InitialCycleNumber'])
                    file.close()
                    if self.verbose:
                        print "wrote profile file ",filename
            del obj
            if callbacks:
                for call in ensure_list(callbacks):
                    profile.add_callback(call)
            if title != None:
                if time == True:
                    title += " %0.2e %s"%(self.ds.parameters['InitialTime']*\
                                        self.UnitMultiplier,self.UnitName)
                self.pc.plots[-1]._axes.set_title(title)
                print title
                profile.add_callback(dave_callback.title(title))
            if self.format == None:
                print self.pc.save(frame_template%(self.returnsetnumber(frame)))
            elif self.format == 'dave':
                n_fields=len(fields)
                dsave(self.pc,frame_template%(self.returnsetnumber(frame)),
                      field_name = "%s_"*n_fields%tuple(fields),
                      ds_list=[self.basename],
                      script_name='uber')
            else:
                print self.pc.save(frame_template%(self.returnsetnumber(frame)), format = self.format)
                #print frame_template%(frame)

    def read_clumps_2(self):
        """second implementation of read clumps, to work with Paper14 needs"""
        if self.clump_prefix == None:
            print "please define clump_prefix; not this should have frame info"
            return
        if self.clumps == None:
            self.clumps = {}
        file = open('stuff_monitor.txt','w')
        file.close()
        t0 = time.time()
        for i in self.frames:
            if self.clumps.has_key(i):
                print "Has Key %d"%i
                continue
            self.clumps[i] = {}
            self.fill(i,get_region=False)
            file = open('stuff_monitor.txt','a')
            file.write('%d reading\n %f'%(i,time.time()-t0))
            file.close()
            to_glob = "%s/%s"%(self.directory,self.clump_prefix%i)
            print to_glob
            all_clumps = glob.glob(to_glob)
            for clump_name in all_clumps:
                print "read", clump_name
                try:
                    clump_ind = int(clump_name[-4:])
                except:
                    print "skip"
                    continue
                self.clumps[i][clump_ind] = clump_stuff.load(clump_name)



    def read_clumps(self):
        """reads self.fields clumps from self.directory + self.clump_prefix"""
        if self.clump_prefix == None:
            print "please define clump_prefix"
            return
        if self.clumps == None:
            self.clumps = {}

        file = open('stuff_monitor.txt','w')
        file.close()
        t0 = time.time()
        for i in self.frames:
            if self.clumps.has_key(i):
                print "Has Key %d"%i
                continue
            self.fill(i,get_region=False)
            #basename = "%s/DD%04d/data%04d"%(self.directory,i,i)
            #exec("ds = %s('%s')"%(self.OutputName,basename))
            ds=self.ds
            print "CURRENT HASH:",ds._hash()
            file = open('stuff_monitor.txt','a')
            file.write('%d reading\n %f'%(i,time.time()-t0))
            file.close()
            filename = "%s/%s_%04d"%(self.directory,self.clump_prefix,i)
            print "=============== %s ==============="%filename
            if glob.glob(filename) == []:
                print "File not found!", filename
                continue
            try:
                self.clumps[i] = clump_stuff.load(filename)
            except KeyError as e:
                print "Key Error:  Probably a bad hash.  Get CTID in the dang code."
                print "You might want to try:"
                print "sed -i 's,%s,%s,g' %s"%(e.args[0][0],ds._hash(),filename)
                raise
            except:
                print "Probably the module error."
                print "sed -i 's,lagos.BaseDataTypes,data_objects.data_containers,g' "
                raise

            self.clumps[i].add_subset(clump_subset.alpha_bound)
            self.clumps[i].add_subset(clump_subset.all_valid)
            self.clumps[i].add_subset(clump_subset.energetically_bound)

    #checkrel def grid(self,n,g,field='Density',gz=0):
    def grid(self,n,g,field='Density',gz=0):
        """Reads dataset *field* from set *n* grid *g*.
        Removes *gz* ghost_zones from either side"""

        #self.set_filename_preset_syntax(n)
        #self.set_filename(n)
        #gridname= self.basename + ".grid%04d"%g
        self.fill(n)
        fname = self.h.grids[g].filename
        fptr = h5py.File(fname,'r')
        try:
            datacube = fptr[field]
        except:
            try:
                gridname = "Grid%08d"%self.h.grids[g].id
                datacube=fptr[gridname][field]
            except:
                raise

        sh=datacube.shape
        dest=np.zeros( sh )
        datacube.read_direct(dest)
        dest = dest.transpose()
        if len(sh) == 1:
            return dest
        else:
            sh=dest.shape
            return dest[gz:sh[0]-gz,gz:sh[1]-gz,gz:sh[2]-gz]
        fptr.close()
    def grid_div_b(self,n,g):
        bx = self.grid(n,g,'BxF')
        #davetools.stat(bx,'bx')
        by = self.grid(n,g,'ByF')
        #davetools.stat(by,'by')
        bz = self.grid(n,g,'BzF')
        #davetools.stat(bz,'bz')
        sh=bx.shape
        dest=np.zeros( sh - nar([1,0,0]))
        dest += bx[1:,:,:]-bx[:-1,:,:]
        dest += by[:,1:,:]-by[:,:-1,:]
        dest += bz[:,:,1:]-bz[:,:,:-1]
        return dest

    def stat(self,field,frame=None):
        """print min and max of *field*"""
        self.fill(frame)
        minTuple = self.h.find_min(field)
        maxTuple = self.h.find_max(field)
        return {'min':minTuple,'max':maxTuple}

    def mstat(self,fields=None,frames=None, format = "%9.2e",Norm=False):
        """Min and max for all fields in self.fields.
        Field list overridden by *fields*.
        *Norm* subtracts off the volume-weighted mean."""
        print fields
        for frame in frames:
            if fields is None:
                fields = self.fields
            self.fill(frame)
            out = self.region.quantities['Extrema'](fields)
            if Norm is True:
                for n, field in enumerate(fields):
                    avg = self.region.quantities['WeightedAverageQuantity'](field,'CellVolume')
                    out[n] = out[n][0]-avg, out[n][1]-avg
            if hasattr(Norm,'has_key'):
                for n, field in enumerate(fields):
                    avg =  0
                    if Norm.has_key(field):
                        avg = Norm[field]
                        print field,avg
                    out[n] = out[n][0]-avg, out[n][1]-avg
            format_string = "%s %s %s"%(format,format,"%s")
            for n, field in enumerate(fields):
                print format_string%(out[n][0], out[n][1], field)

    def vstat(self,field,frame=None, grid_list=None):
        """Prints stats on each grid.  For the impatient."""
        self.fill(frame)
        if grid_list is None:
            grid_list = slice(None)
        for i,g in enumerate(self.h.grids[grid_list]):
            this_max = g[field].max()
            this_min = g[field].min()
            if i == 0:
                minimum = this_min
                maximum = this_max
            minimum = min(this_min, minimum)
            maximum = max(this_max, maximum)
            print "%4d, (%0.2e, %0.2e) (%0.2e, %0.2e)"%(g,this_min,this_max,minimum,maximum)


    def save_set(self,set,frame=None,field=None,prefix="cube", num_ghost_zones=0, debug=0, dtype='float32'):
        """Saves *set* to disk using standard naming conventions."""
        if frame == None:
            frame = self.frames[0]
        if field == None:
            field = self.fields[0]
        directory = self.product_dir(frame)
        filename= "%s/%s_%s.%s.%d"%(directory,prefix,field,dtype,num_ghost_zones)
        if debug > 0:
            print "saving", filename
        file = h5py.File(filename,'w')
        file.create_dataset(field,set.shape, data=set, dtype=dtype)
        file.close()




        
