import yt
import numpy as np
import pdb
import pylab
import matplotlib.pyplot as plt
import matplotlib
import copy
nar = np.array
def stat(array,strin='', format='%0.16e'):
    template = '['+format+','+format+'] %s %s'
    print(template%(array.min(),array.max(),array.shape,strin))

def dumb_map(V,cmap = 'jet',range=[0,1]):
    """Does all the crap that matplot lib needs to map V, assumed to be in [0,1] to the colormap"""
    norm = matplotlib.colors.Normalize(range[0],range[1])
    exec('cmap = matplotlib.cm.%s'%cmap)
    color_map = matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap)
    return color_map.to_rgba(V)
        
def tube(ds_list,fields=None,times=None, points=[(0.0,0.0,0.0),(1.0,0.0,0.0)],width=2,
         filename='FILE.png',fig_args={},plot_args={},legend=False,
         renorm=False,return_ray=False,
         offsets=None,
         axis=0,coord=(0.505,0.505), coord2=None,ray_type='oray',
         debug=0,delta=False,
         labels=None, xlim = None, ylim={}, ylim_method=None, yscale={}): 
    """Plots 1d profiles for each ds in *ds_list*, for each field in *fields* and time in *times*.
    One plot per field, with oober and time overplotted in each field window.
    Uses the AMRRay from yt, with start and end points given by *points*.
    If *fields* is not given, the minimum set is taken from oober.fields for each oober in *ooberlist*.
    If *times* is not given, each oober uses the values in oober.frames
    *renorm* removes the mean of each field.
    *offsets*, if not None, should be a list of multiplicitive offsets for each oober
    *axis*,*coord* gives the arguments to ortho_ray
    *delta* = True takes the difference of oober[0] with oober[1:]


    usage:
    frame_list = [0,1,2,3]
    sim_5="/whereIputit"
    ds_list = [yt.load('%s/DD%04d/data%04d'%(sim_5,n1,n1)) for n1 in frame_list]
    fields = ['density','pressure']
    renorm={'density':1,'pressure':0.6,'TotalEnergy':0.9}
    y=tube.tube(ds_list,  delta=True, fields=fields, renorm=renorm,filename = "thing.png")
    """
    #global plt
    #ds_list = ensure_list(ds_list)
    fig = plt.figure()
    ray_set = {}

    # Set up times and fields to loop over

    if fields == None:
        fields = ['density']
    if legend == 1: 
        extra_for_legend = 1
    else:
        extra_for_legend = 0
    n_rows = int(1.*(len(fields)+extra_for_legend)/width+0.5)
    # Make list of rays
    #pdb.set_trace()
    if 0:
        """this is a dirty hack."""
        if labels is None:
            labs = list(range(len(ds_list)))
        else:
            labs = labels
        for ds, lab in zip(ds_list, labs):
            ray_set[lab]= ds.ortho_ray(axis,coord) 
        return ray_set
    else:
        xcoord={}
        for nds,ds in enumerate(ds_list):
            if ray_type == 'oray':
                ray_set[ds]= ds.ortho_ray(axis,coord) 
            elif ray_type == 'obl':
                if len(coord) != 3 or len(coord2) != 3:
                    print("coordinates must be 3d")
                    return -1
                ray_set[ds] = ds.r[ coord:coord2 ]
            elif ray_type == 'multi':
                ray_set[ds] = ds.r[ coord[nds]:coord2[nds]]
            if type(axis) == str: 
                xcoord[ds]=axis
            elif type(axis) == int:
                xcoord[ds]='xyz'[axis]    
            elif type(axis) == list:
                xcoord[ds] = axis[nds]

    #make the figure
    fig = plt.figure(**fig_args)
    first_ax = None
    Handels=[]
    
    units={'density':'code_mass/code_length**3','pressure':'code_mass/(code_length*code_time**2)'}
    for x in 'xyz':
        units['B'+x] = 'code_magnetic'
        units[x+'-velocity']='code_length/code_time'

    for i,field in enumerate(fields):

        fiducial_field = []
        ax = fig.add_subplot(n_rows,width,i+1)
        if debug > 0:
            print("n_rows %d width %d i %d field %s"%(n_rows,width,i+1, field))
        counter = -1
        for n_ds,ds in enumerate(ds_list):
            counter += 1
            this_ray = ray_set[ds]
            this_x = copy.copy(this_ray[xcoord[ds]].v[:])
            sort_x = np.argsort(this_x)
            this_x = this_x[sort_x]
            this_y = copy.copy(this_ray[field])
            if field in units:
                this_y = this_y.in_units(units[field])
            this_y =this_y.v[sort_x]
            #pdb.set_trace()
            if delta:
                if n_ds==0:
                    fiducial_field = this_y
                    continue
                else:
                    this_y -= fiducial_field
                    nonzero = fiducial_field != 0
                    this_y[nonzero] /= fiducial_field[nonzero]
            stat(this_y,"before %s %s %0.2f"%(ds,field,ds['InitialTime']), format="%0.2e")
#           if renorm is True:
#               this_y = this_y - (this_y*this_x).sum()/this_x.sum()


            if hasattr(renorm,'has_key'):
                if field in renorm:
                    this_y -= renorm[field]
            if offsets is not None:
                this_y *= offsets[n_ds]
            #stat(this_y,"after %s %s t=%0.2f"%(field,labels[n_ds],ds['InitialTime']),format='%0.2e')


            #plot_args['color'] = dumb_map(counter, range=[0,len(time_set[o.name])])
            if first_ax is None:
                plot_args['label'] =  '%s %0.2f'%(labels[n_ds],ds['InitialTime'])
            #plot_args['linewidth']=0.1
            #print this_y
            L = ax.plot(this_x,this_y,**plot_args)
            if ylim_method == 'monotonic':
                if field not in ylim:
                    ylim[field] = [min(this_y),max(this_y)]
                else:
                    ylim[field][0] = min( ylim[field][0], min(this_y))
                    ylim[field][1] = max( ylim[field][1], max(this_y))
                print(("YLIM",ylim))

            if this_y.min() > 0 :
                ax.set_yscale( yscale.get(field,'log'))

#               if first_ax is None:
#                   Labels.append( '%s %d'%(o.name,t))

        if first_ax is None:
            first_ax = ax
            Handels, LabelsToUse = first_ax.get_legend_handles_labels()
        if xlim is not None:
            ax.set_xlim(xlim)
        if field in ylim:
            ax.set_ylim(ylim[field])

        ax.set_ylabel(field)

    if legend == 1:
        ax = fig.add_subplot(n_rows,width,i+2)
        ax.legend( Handels, LabelsToUse,loc = 0)
    elif legend == 2:
        first_ax.legend(Handels,LabelsToUse)
    fig.savefig(filename)
    plt.close(fig)
    print(filename)
    if return_ray:
        print("There's something persistent that's getting messed up with returning the ray. Be careful.")
        return ray_set
    else:
        return None
