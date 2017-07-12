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
nar = np.array
array=np.array
import pdb
import os
import h5py
import types, time,weakref,davetools
#import dave_callbacks
from davetools import dsave, no_trailing_comments, ensure_list
import time, fPickle, glob
if not_ported:
    import clump_stuff, clump_subset
    from clump_list_sort import *
    #reload(clump_stuff)
import matplotlib.colorbar as cb
import re
import copy

class fleet():
    def __init__(self,taxi_list=[]):
        self.taxi_list = []
        self.next_taxi_index=0
        for car in taxi_list:
            if isinstance(car,types.StringType):
                self.taxi_list.append(taxi(car))
            else:
                self.taxi_list.append(car)
        self.namelength = max([len(t.name) for t in self.taxi_list])
        self.nt = "%"+str(self.namelength+1)+"s "

    def next(self):
        this_taxi_index = self.next_taxi_index
        self.next_taxi_index += 1
        if self.next_taxi_index >= len(self.taxi_list) + 1:raise StopIteration
        return self.taxi_list[this_taxi_index]
    def __iter__(self):
        self.next_taxi_index=0
        return self
    def __getitem__(self,item):
        out = []
        if isinstance(item,types.IntType):
            out = self.taxi_list[item]
        else:
            for car in self.taxi_list:
                print self.nt%car.name, car.__dict__[item]
                out.append(car.__dict__[item])
        return out
    def __setitem__(self,item,value ):
        for car in self.taxi_list:
            car.__dict__[item] = value
            
    def outname(self,prefix):
        """sets car.outname = prefix+car.name"""
        for car in self.taxi_list:
            car.outname = prefix + car.name
    def plot(self,*args, **kwargs):
        #check for fixed colorbars.
        prior_zlim=None
        if kwargs.has_key('prefix'):
            self.outname(kwargs.pop('prefix'))

        for car in self.taxi_list:
            if car.Colorbar is 'fixed' and prior_zlim is not None and car.zlim is not None:
                car.zlim = prior_zlim
            car.plot(*args, **kwargs)
            prior_zlim = car.zlim
    def __call__(self,string, frames=False):
        """Execute arbitrary code on the cars in the taxi fleet.
        output can be used to return values."""
        output = []
        for car in self.taxi_list:
            if frames is False:
                exec(string)
            else:
                for frame in frames:
                    car.fill(frame)
                    exec(string)
        return output
    def save(self,suffix=""):
        for car in self.taxi_list:
            thisname = car.name+suffix
            car.save(thisname)
    def allnames(self):
        ncar=len(self.taxi_list)
        return "%s_"*ncar%tuple([car.name for car in self.taxi_list])
    def dumb_profile(self,fields,**kwargs):
        for car in self.taxi_list:
            car.dumb_profile(fields,**kwargs)
    def phase(self,*args,**kwargs):
        for car in self.taxi_list:
            car.phase(*args,**kwargs)
        
    def profile(self,*args,**kwargs):
        """Runs profile on all using *args and **kwargs.
        Also plots the combined set"""
        for car in self.taxi_list:
            car.profile(*args,**kwargs)
        plt.clf()
        color_by = kwargs.pop('color_by','sim') #either color with sim or frame. The other gets line style
        ntotal = max([len(car.frames) for car in self.taxi_list])
        simtotal = len(self.taxi_list)
        if color_by=='frame':
            rm = davetools.rainbow_map(ntotal)
        elif color_by=='sim':
            rm = davetools.rainbow_map(simtotal)
        plt.clf()
        linelist = {0:'-',1:'--',2:':'}
        all_frames=[]
        plot_args={'linewidth':0.3}  #I can probably make this the same as kwargs
        prior_y = None
        for n,car in enumerate(self.taxi_list):
            all_xbins = car.profile_data['all_xbins']
            all_profiles = car.profile_data['all_profiles']
            for i,frame in enumerate(car.frames):
                thex= all_xbins[i]
                they=all_profiles[i]
                plot_args['label']='%s %04d'%(car.name,frame)
                if prior_y is not None and False:
                    l1norm = np.mean(np.abs( they - prior_y ))/np.mean(they)
                    print l1norm
                    prior_y = they
                    if l1norm > 1e-8:
                        they *= 1.1
                        plot_args['label'] += "+off"
                else:
                    prior_y = they
    
                if color_by=='sim':
                    plot_args['c']=rm(n)
                    plot_args['linestyle']=linelist.get(i,'-')
                elif color_by=='frame':
                    plot_args['c']=rm(i)
                    plot_args['linestyle']=linelist.get(n,'-')
                plt.plot( thex, they,**plot_args)
                #plt.plot( all_xbins[i], all_profiles[i],c=rm(i),linestyle=linelist[n],label='%s %04d'%(car.name,frame))
                if frame not in all_frames:
                    all_frames.append(frame)
        plt.xlabel(car.profile_data['fields'][0])
        plt.ylabel(car.profile_data['fields'][1])
        plt.xscale(car.profile_data['scales'][0]); plt.yscale(car.profile_data['scales'][1])
        plt.legend(loc=0)
        frame_name = "_%04d"*len(all_frames)%tuple(all_frames)
        profname = '%s_prof_%s_%s_n%s.pdf'%(self.allnames(), car.profile_data['fields'][0], car.profile_data['fields'][1], frame_name)
        print profname
        plt.savefig(profname)

    def find_extrema(self,fields=None,frames=None, manual_positive=False):
        all_fields = []
        for car in self: #I feel like there's an easier way to do this.
            for field in car.fields:
                if field not in all_fields:
                    all_fields.append(field)

        extrema_store={}
        for n,car in enumerate(self.taxi_list):
            car.find_extrema(fields,frames, manual_positive=manual_positive)
            for field in all_fields:
                if car.extrema.has_key(field):
                    if extrema_store.has_key(field):
                        extrema_store[field][0] = min([extrema_store[field][0],car.extrema[field][0]])
                        extrema_store[field][1] = max([extrema_store[field][1],car.extrema[field][1]])
                    else:
                        extrema_store[field]=np.zeros(2)
                        extrema_store[field][0] = car.extrema[field][0]
                        extrema_store[field][1] = car.extrema[field][1]
        for car in self.taxi_list:
            car.extrema = copy.copy(extrema_store)







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
    taxi.read(file) populates itself from file by executing each line in the file.
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

        #List of members that will NOT be written.  
        self.ExcludeFromWrite = ['ExcludeFromWrite']
        

        #To make the file easier to use, some of the uber members are written first
        self.WriteMeFirst = ['name','directory','outname','operation','frames','fields','axis',
                            'name_syntax', 'name_files','name_dir','GlobalParameters',
                             'ProfileDir','ProfileName']
        self.ExcludeFromWrite.append('WriteMeFirst')
        self.ExcludeFromWrite.append('profile_data')
        self.ExcludeFromWrite.append('last_prof')

        #Some members may need special treatment.
        self.WriteSpecial = {}
        self.ExcludeFromWrite.append('WriteSpecial')
        #externals
        self.OutputLog = "OutputLog"
        self.plot_origin = 'domain'  #see yt docs.
        self.axis = 0                #Axis for the plot
        self.frames = []             #Dump numbers to be plotted.
        self.fields = 'Density'      #Field (hopefully to be list) for the plot
        self.weight_field = None           #Weight field.  Projections?    
        self.name   = 'sim'          #identifier for this taxi instance
        self.outname = 'Image'       #Prefix for output
        self.directory = '.'         #Directory where the simulation was run. (Where the DD* dirs are)
        self.ProfileDir= "./ProfileFiles/"
        self.ProfileName = None
        self.set_log = {}
        self.subdir = True           #Are there DD0000 directories, or is it just data0000.grid*?
        self.operation = 'Full'      #The operation that will be done.  taxi.operations() for a list.
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
        self.last_particle_index_step = None
        self.ExcludeFromWrite.append('last_particle_index_step')
        self.last_particle_indices = None
        self.ExcludeFromWrite.append('last_particle_indices')
        self.field_parameters = {}
        self.left=np.array([0.0,0.0,0.0])      #The left of the region extraction.
        self.right=np.array([1.0,1.0,1.0])     #The right of the extracted region.
        self.radius=None
        self.height=None
        self.restrict = False                   #Restric plot to img_left, img_right
        self.img_left = None                   #left edge of image
        self.img_right = None                  #right edge of image.
        self.img_width = None                  #width of image.
        self.width = None
        self.zoom_sequence=None
        self.region_type = 'all'
        self.center= 0.5*(self.left+self.right)#Center, used for all things that need a center
        self.periodic = False              #For projection shift
        self.cmap = {'default':None}       #color map, field driven. Everyone uses default.
        self.Colorbar = None               #Monotone: monotonically increaseing, forced to min,max
        self.use_colorbar = True           #use color bar for plots.  Bad naming.
        self.zlim = {}
        self.slice_zlim = {}
        self.proj_zlim = {}
        self.callbacks = []                #Which callbacks to use.  taxi.callbacks() for options
        self.callback_args={'particles':{'args':[1],'kwargs':{'col':'r'}},'select_particles':{'args':1,'kwargs':{'col':'y'}}}
        self.callback_args['star_particles']={'args':[1],'kwargs':{'col':'g'}}
        self.callback_args['nparticles']= {'args': [[0.03, 0.03]], 'kwargs': {'coord_system': 'axis', 'text_args': {'color': 'red'}}}
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
        self.name_syntax = 'outputlog'  #outputlog or preset
        self.name_files = "data"        #prefix on data file
        self.name_dir = "DD"  #Prefix on said directories
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
        self.ExcludeFromWrite.append('the_plot')
        self.hold_values={}

        self.particle_filter_names=[]  #Filters get added directly to the ds, e.g.
        #http://yt-project.org/doc/cookbook/calculating_information.html?highlight=formed_star

        self.derived_fields={} #Derived fields can also be added to the ds directly.
        #This is not fully developed, syntax should be sth like
        #self.derived_fields['field_name'] = {'function':{dict of args}}
        #which then gets called on ds as its loaded.
 


        if filename != None:
            #If no directories are mentioned, assume the actual file is ./taxi_stand/filename
            try:
                filename.index('/')
                actual_filename = filename
            except:
                actual_filename = "taxi_stand/%s"%filename
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
        try:
            filename.index('/')
            actual_filename = filename
        except:
            actual_filename = "taxi_stand/%s"%filename
        file = open(actual_filename,"w")
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
        print "saved", actual_filename

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
            frame = self.return_frames()[-1]
        #take 1: no subdirectory (very old style)
        self.basename_template = self.directory+'/'+self.name_files+'%04i'
        self.basename = self.basename_template %(frame)
        if glob.glob(self.basename+".hierarchy") == []:
            self.basename_template = self.directory + "/" + self.name_dir + '%04i/'\
                                     + self.name_files + "%04i"
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
                frame = self.return_frames()[0]
            return self.frame_dict[frame]['SetNumber']

    def get_frames(self):
        """parses OutputLog to populate self.frame_dict"""
        DirSetNumber = re.compile(r'([^\d]*)(\d\d\d\d)/([^\d]*)(\d\d\d\d)')
        SetNumber = re.compile(r'([^\d]*)(\d\d\d\d)')
        SetNumber = re.compile(r'(.*)(\d\d\d\d)')
        debug = 0
        if not hasattr(self,'frame_dict') or self.frame_dict is None:
            self.frame_dict = {}
            output_name = "%s/%s"%(self.directory , self.OutputLog)
            output_log = open(output_name, 'r')
            print "OutputName", self.name, output_name

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
                    self.frame_dict[nframe]['name_files'] = match.group(3)
                    self.frame_dict[nframe]['SetNumber'] = int(match.group(2))
                elif match2 is not None:
                    self.frame_dict[nframe]['DirPrefix']=None
                    self.frame_dict[nframe]['name_files'] = match2.group(1)
                    self.frame_dict[nframe]['SetNumber'] = int(match2.group(2))
                else:
                    print "error in set parsing"
                    print linesplit[2]

                    self.frame_dict[nframe]['DirPrefix']=None
                    self.frame_dict[nframe]['name_files']  =None
                    self.frame_dict[nframe]['SetNumber']=-1

                try:
                    self.frame_dict[ nframe ]['time']   = float(linesplit[4])
                except:
                    self.frame_dict[ nframe ]['time']   = -1
                #pdb.set_trace()
    def set_filename_from_outputlog(self,frame=None):
        self.get_frames()
        if frame == None:
            frame = self.return_frames()[-1]
        self.basename = "%s/%s"%(self.directory,self.frame_dict[frame]['dsname'])
        if self.verbose == True:
            print "==== %s ===="%(self.basename)

    def fill(self, frame = None):
        """populate parameter file, plot collection, hierarchy, region if desired.
        Possibly could be renamed 'load' """
        self.set_filename(frame)
        self.ds = yt.load(self.basename)
        for filter_name in self.particle_filter_names:
            self.ds.add_particle_filter(filter_name)
        #na_errors= np.seterr(all='ignore')
        #np.seterr(**na_errors)

    def get_region(self, frame=None):
        if frame is not None or self.ds is None:
            self.fill(frame)
        if self.region_type.lower()=='sphere':
            reg = self.ds.sphere(self.center, self.radius)
        if self.region_type.lower() in ['rectangle','region']:
            self.center = 0.5*(self.left + self.right)
            reg = self.ds.region(self.center, self.left,self.right)
        if self.region_type.lower() in ['all','all_data']:
            reg = self.ds.all_data()
        if self.region_type.lower() in ['disk']:
            print "OK!"
            reg = self.ds.disk(self.center,self.normal,self.radius,self.height)
        if self.region_type.startswith('grid'):
            """gridN gets grid N."""
            N = int( self.region_type[4:])
            g = self.ds.index.grids[N]
            center = 0.5*(g.LeftEdge+g.RightEdge)
            reg = self.ds.region(center,g.LeftEdge,g.RightEdge)
        for parameter in self.field_parameters.keys():
            reg.set_field_parameter(parameter, self.field_parameters[parameter])
        #self.reg  = weakref.proxy(reg)
        return reg




    def arg_setter(self,kwargs):
        for arg in kwargs.keys():
            if arg.startswith("_"):
                continue
            if self.__dict__.has_key(arg):
                self.__dict__[arg] = kwargs[arg]


    def return_frames(self):
        if not hasattr(self,'frame_dict') or self.frame_dict is None:
            self.get_frames()
        all_frames = sorted(self.frame_dict.keys())
        if self.frames=='all':
            return_frames = all_frames
        elif self.frames.startswith('every'):
            interval = int(self.frames.split(" ")[-1])
            return_frames = all_frames[::10]
            if return_frames[-1] != all_frames[-1]:
                return_frames += all_frames[-1:]

        elif self.frames == 'last':
            return_frames = all_frames[-1]
        elif self.frames == 'all_reverse':
            return_frames = all_frames[::-1]
        else:
            return_frames = self.frames

        return return_frames

    def plot(self,local_frames=None, **kwargs):
        self.arg_setter(kwargs)
        if 'new_particles' in self.callbacks:
            if len(self.axis) *len(self.fields) > 1:
                print "WARNING: new_particles callbacks does not work with multiple axes or fields."
        print self.zlim
        """Makes the uber plot.  Options for plot can be found in the uber description,
        plot operations can be found by typing uber.operations()"""
        start_time = time.time()    
        
        frame_template = self.outname + "_%04i"
        FirstOne = True #For monotonic plots.

        if local_frames==None:
            local_frames=self.return_frames()

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

                    if self.set_log.has_key(field):
                        plot_or_fig.set_log( field, self.set_log[field])

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
                        plot_or_fig.set_origin('domain')
                        if self.zoom_sequence is not None:
                            for nz,width in enumerate(self.zoom_sequence):
                                plot_or_fig.set_width(width)
                                this_frame_template = frame_template + "_zoom%02d"%nz
                                print plot_or_fig.save(this_frame_template%frame)
                        else:
                            print plot_or_fig.save(frame_template%frame)
    def set_center_max(self,field='density', frame=None):
        if self.ds is None or frame is not None:
            self.fill(frame)
        value,center = self.ds.find_max(field)
        self.center = np.array(center)
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
                the_plot = yt.ProjectionPlot(self.ds,axis,field, center=self.center, weight_field=self.weight_field, **extra_args)
                #self.pc.add_projection(field,axis, use_colorbar=use_colorbar,
                #                            **extra_args)
                self.zlim = self.proj_zlim

            elif self.operation == 'RegionProjection_kill_this':
                """Might be broken"""
                #setup a projection
                self.center = 0.5*(self.left + self.right)
                reg = self.region(self.center, self.left,self.right,field)
                self.reg = weakref.proxy(reg)
                self.proj = self.ds.proj(axis,field,center=self.center,source=self.reg,periodic=self.periodic,weight_field=self.weight_field)
                the_plot = self.proj.to_pw()
                self.zlim = self.proj_zlim
            elif self.operation == 'RegionProjection':
                #setup a projection
                reg = self.get_region()
                self.proj = self.ds.proj(field,axis,center=self.center,data_source=reg, weight_field=self.weight_field) #,periodic=self.periodic)
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
                do_log = True #please get this from the yt object.
                ok_zones = np.isnan(the_plot.data_source[field]) == False
                if do_log:
                    ok_zones = np.logical_and(ok_zones, the_plot.data_source[field] != 0)
                if ok_zones.sum() > 0:
                    this_min = the_plot.data_source[field][ok_zones].min()
                    this_max = the_plot.data_source[field][ok_zones].max()
                if self.Colorbar in ['Monotone', 'monotonic']:
                    if FirstOne:
                        self.zlim[field] = [this_min,this_max]
                    else:
                        self.zlim[field][0] = min([self.zlim[field][0],this_min])
                        self.zlim[field][1] = max([self.zlim[field][1],this_max])
                    if self.verbose:
                        print "set lim", self.zlim[field]
                elif self.Colorbar in  ['Fixed', 'fixed'] :
                    if not self.zlim.has_key(field):
                        do_log = True #please get this from the yt object.
                        self.zlim[field] = [this_min,this_max]
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
            raise
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
                                                   self.name_dir,
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
            frame = self.return_frames()[0]
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
            frame = self.return_frames()[0]
        if field == None:
            field = self.fields[0]
        directory = self.product_dir(frame) #"%s/%s%04d.products/"%(self.directory,self.name_dir,frame)
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
    def find_extrema(self,fields=None,frames=None,manual_positive=False):
        if fields is not None:
            local_fields = fields
        else:
            local_fields = self.fields
        if frames:
            local_frames = frames
        else:
            local_frames = self.return_frames()
        for frame in local_frames:
            self.fill(frame)
            reg=self.get_region(frame)
            for field in local_fields:
                if manual_positive:
                    vals = reg[field]
                    positive = vals > 0
                    this_extrema=[0,0]
                    this_extrema[0] = vals[positive].min()
                    this_extrema[1] = vals[positive].max()
                else:
                    this_extrema = reg.quantities['Extrema'](field,non_zero=True)
                if not self.extrema.has_key(field):
                    self.extrema[field] = [this_extrema[0].v, this_extrema[1].v]
                else:
                    self.extrema[field][0] = min([this_extrema[0].v, self.extrema[field][0]])
                    self.extrema[field][1] = max([this_extrema[1].v, self.extrema[field][1]])

    def dumb_profile(self,fields,**kwargs):
        frame_template = self.outname + "_%04i"
        for frame in self.return_frames():
            reg = self.get_region(frame)
            plot=yt.ProfilePlot(reg,fields[0],fields[1],**kwargs)
            plot.save(frame_template%frame)
    def profile(self,fields, callbacks=None, weight_field=None, accumulation=False, fractional=True, scales=['log','log'],n_bins=64,extrema=None, units=[None,None]):
        """needs to be generalized with bins."""
        frame_template = self.outname + "_%04i"
        weight_name = ""
        if weight_field is not None:
            weight_name = "%s_"%weight_field
            frame_template += "_%s"%weight_field
        all_xbins = []
        all_profiles = []
        for frame in self.return_frames():
            plt.clf()
            reg = self.get_region(frame)
            local_extrema = None
            prof = yt.create_profile(reg,fields[0],fields[1] ,weight_field=weight_field,accumulation=accumulation,
                                    fractional=fractional, n_bins=n_bins, extrema=extrema)
            self.last_prof=prof
            the_x = 0.5*(prof.x_bins[1:]+prof.x_bins[0:-1])
            the_y = prof[fields[1]]
            if units[0] is not None:
                the_x = the_x.in_units(units[0])
            if units[1] is not None and fractional is not True:
                the_y = the_y.in_units(units[1])
            x_units=the_x.units
            y_units=the_y.units
            plt.plot(the_x,the_y,label="n%04d"%frame)
            all_xbins.append(the_x)
            all_profiles.append(the_y)
            scaledict={True:'log',False:'linear'}
            plt.xscale(scaledict[scales[0]]); plt.yscale(scaledict[scales[1]])
            plt.xlabel(r'%s $%s$'%(fields[0],x_units)); plt.ylabel(r'%s $%s$'%(fields[1],y_units))
            profname = '%s_prof_%s_%s_n%04d.pdf'%(self.outname, fields[0], fields[1], frame)
            plt.legend(loc=0)
            plt.savefig(profname)
            print profname
            if True:
                if glob.glob(self.ProfileDir) == []:
                    os.mkdir(self.ProfileDir)
                    print "made directory", self.ProfileDir
                filename = "%s/%s_%s_%s"%(self.ProfileDir,frame_template%(self.returnsetnumber(frame)),fields[0],fields[1])
                all_files = glob.glob(filename+"*")
                filename += "_%d.h5"%(len(all_files))
                fptr=h5py.File(filename,"w")
                fptr.create_dataset(fields[0],prof.x_bins.shape,data=prof.x_bins)
                fptr.create_dataset(fields[1],prof[fields[1]].shape,data=prof[fields[1]])
                fptr.create_dataset('pf_name',data=self.name)
                fptr.create_dataset('InitialTime',data=self.ds['InitialTime'])
                fptr.create_dataset('InitialCycle',data=self.ds['InitialCycleNumber'])
                fptr.close()
                if self.verbose:
                    print "wrote profile file ",filename
        ntotal = len(self.return_frames())
        rm = davetools.rainbow_map(ntotal)
        plt.clf()
        for i,n in enumerate(self.return_frames()):
            plt.plot( all_xbins[i], all_profiles[i],c=rm(i),label="n%04d"%n)
        scaledict={True:'log',False:'linear'}
        plt.xscale(scaledict[scales[0]]); plt.yscale(scaledict[scales[1]])
        plt.xlabel(r'%s $%s$'%(fields[0],x_units)); plt.ylabel(r'%s $%s$'%(fields[1],y_units))
        plt.legend(loc=0)
        allframes = "_%04d"*ntotal%tuple(self.return_frames())
        profname = '%s_prof_%s_%s_%sn%s.pdf'%(self.outname, fields[0], fields[1],weight_name, allframes)
        print profname
        plt.savefig(profname)
        self.profile_data={'all_xbins':all_xbins,'all_profiles':all_profiles, 'scales':scales, 'fields':fields}

    def phase(self,fields, callbacks=None, weight_field=None, phase_args={},save=True, n_bins=[64,64], phase_callbacks=[]):
        """Uber wrapper for phase objects.
        for all frames in self.frame, run a phase plot object on the region.
        Run all callbacks on the plot.
        *save* pickles the region in PhaseFiles"""
        if len(fields) != 3:
            print "wrong number of fields: needs 3.", fields
            return
        frame_template = self.outname + "_%04i"
        phase_list = []
        for frame in ensure_list(self.return_frames()):
            reg = self.get_region(frame)
            local_extrema = None
            if self.Colorbar in ['fixed', 'monotonic']:
                if self.extrema.has_key(fields[0]) and self.extrema.has_key(fields[1]):
                    local_extrema = {fields[0]:self.extrema[fields[0]], fields[1]:self.extrema[fields[1]]}
                else:
                    local_extrema = None

            print "LOCAL", local_extrema
            phase_args['bin_fields']=[fields[0],fields[1]]
            phase_args['fields']=[fields[2]]
            phase_args['weight_field']=weight_field
            phase_args['extrema']=local_extrema
            phase_args['n_bins']=n_bins
            phase = yt.create_profile(reg,**phase_args)
            #self.phase = weakref.proxy(phase)
            pp = yt.PhasePlot.from_profile(phase)
            pp.set_xlabel(fields[0])
            pp.set_ylabel(fields[1])
            #print pp.save('derp3.png')
            outname = "%s_%04d"%(self.outname,frame)

            if self.Colorbar in ['fixed','monotonic']:
                xmin=phase.x_bins[0]
                xmax=phase.x_bins[-1]
                ymin=phase.y_bins[0]
                ymax=phase.y_bins[-1]
                if self.extrema.has_key(fields[0]) and self.Colorbar == 'monotonic':
                    xmin=min([self.extrema[fields[0]][0], phase.x_bins[0]])
                    xmax=max([self.extrema[fields[0]][1], phase.x_bins[-1]])
                if self.extrema.has_key(fields[1]) and self.Colorbar == 'monotonic':
                    ymin=min([self.extrema[fields[1]][0], phase.y_bins[0]])
                    ymax=max([self.extrema[fields[1]][1], phase.y_bins[-1]])
                if self.Colorbar in ['fixed', 'monotonic']:
                    if not self.extrema.has_key(fields[0]):
                        self.extrema[fields[0]] = [xmin,xmax]
                    if not self.extrema.has_key(fields[1]):
                        self.extrema[fields[1]] = [ymin,ymax]
                print "extrema, post", self.extrema
                pp.set_xlim( lim_down(self.extrema[fields[0]][0]), lim_up(self.extrema[fields[0]][1]))
                pp.set_ylim( lim_down(self.extrema[fields[1]][0]), lim_up(self.extrema[fields[1]][1]))

            for callback in phase_callbacks:
                callback(pp)
                #key = pp.plots.keys()[0]
                #this_axes = pp.plots[key].axes
                #pp.save('horm.png')
                ##this_axes.plot([1e7,1e10],[1e7,1e10])
                #this_axes.plot([1e-25,1e-23],[1e-25,1e-23])
                #pp.save('derp.png')

            phase_list.append(pp)
            print pp.save(outname)
            #
            # some old hdf5 and callback code. Harvest later.
            #
            if 0:
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
            if 0:
                for call in ensure_list(callbacks):
                    phase.add_callback(call)

        return phase_list
    def count_particles(self, frame=None, ptype='all'):
        """Turns out Metadata.NumberOfParticles isn't always updated close to the output in Enzo."""
        """Also, if ds is loaded already, just leave frame as None and it won't reload"""
        nparticles = 0
        reg = self.get_region(frame)
        try: 
            nparticles = reg[ptype,'particle_index'].size
        except:
            pass
        return nparticles

    def get_new_indices(self):
        reg=self.get_region()
        try:
            these_indices = reg['particle_index']
        except:
            these_indices = self.ds.arr([],'dimensionless')
        if self.last_particle_indices is not None and self.last_particle_index_step is not self.ds['InitialCycleNumber']:
            to_plot = np.setdiff1d(these_indices,self.last_particle_indices)
        else:
            to_plot=these_indices
        self.last_particle_indices=these_indices
        self.last_particle_index_step = self.ds['InitialCycleNumber']
        return to_plot
    def add_callbacks(self,the_plot,axis=None):
        for i,callback in enumerate(self.callbacks):
            args = self.callback_args.get(callback,{'args':[],'kwargs':{}})
            myargs=args.get('args',[])
            mykwargs=args.get('kwargs',{})
            if isinstance(callback,types.StringType):
                if callback == 'velocity':
                    the_plot.annotate_velocity()
                elif callback == 'grids':
                    the_plot.annotate_grids()
                elif callback == 'particles':
                    nparticles = self.count_particles()
                    if nparticles>0:
                        the_plot.annotate_particles(*myargs,**mykwargs)
                elif callback == 'text':
                        the_plot.annotate_text(myargs[0],myargs[1],**mykwargs)
                elif callback == 'nparticles':
                        the_plot.annotate_text(myargs[0],r'$n_p=%d$'%self.count_particles(),**mykwargs)
                elif callback == 'star_particles':
                    nparticles = self.count_particles()
                    if nparticles>0:
                        reg = self.get_region()
                        these_indices = reg['formed_star','particle_index']
                        mykwargs['ptype']='formed_star'
                        the_plot.annotate_select_particles(myargs[0], **mykwargs)
                        pargs = self.callback_args['nparticles']['args']
                        pkwargs=self.callback_args['nparticles']['kwargs']
                        the_plot.annotate_text(pargs[0],r'$n_{\rm{new}}=%d$'%these_indices.shape,**pkwargs)
                elif callback == 'time_title':
                    units  = self.callback_args.get('time_title',{'units':'Myr'}).get('units','Myr')
                    format = self.callback_args.get('time_title',{'format':'%0.2f'}).get('format','%0.2f')
                    output = r"$"
                    time = self.ds.current_time.in_units(units)
                    output += format%time
                    output += "\ %s$"%time.units
                    the_plot.annotate_title(output)
                elif callback == 'new_particles':
                    """This callback does not work with multiple plots for each frame.
                    This is due to the fact that the old particle list is updated each call
                    A more sophisiticated cache system, or perhaps updating the 'last particle list'
                    in the plot loop is necessary."""
                    nparticles = self.count_particles()
                    if nparticles>0:
                        these_indices = self.get_new_indices()
                        mykwargs['indices'] = these_indices
                        mykwargs['col'] = 'r'
                        the_plot.annotate_select_particles(myargs[0], **mykwargs)
                        pargs = self.callback_args['nparticles']['args']
                        pkwargs=self.callback_args['nparticles']['kwargs']
                        the_plot.annotate_text(pargs[0],r'$n_{\rm{new}}=%d$'%these_indices.shape,**pkwargs)

                elif callback == 'magnetic_field':
                    the_plot.annotate_magnetic_field()
                else:
                    raise PortError("Callback %s not supported"%callback)
            else:
                print "Where the heck did you get that callback at?"
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
            for n in self.return_frames():
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
            frames = self.return_frames()

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
        for frame in ensure_list(self.return_frames()):
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
        for i in self.return_frames():
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
        for i in self.return_frames():
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
            frame = self.return_frames()[0]
        if field == None:
            field = self.fields[0]
        directory = self.product_dir(frame)
        filename= "%s/%s_%s.%s.%d"%(directory,prefix,field,dtype,num_ghost_zones)
        if debug > 0:
            print "saving", filename
        file = h5py.File(filename,'w')
        file.create_dataset(field,set.shape, data=set, dtype=dtype)
        file.close()




        
