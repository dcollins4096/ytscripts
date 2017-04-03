#execfile('go_lite')
execfile('go')
import taxi
from p19b_select_particle_callback import *
import time
import clump_properties
reload(clump_properties)
#Clump pre-image analysis
#1.) Run amazing huge sims.  The best.
#2.) Run clumps on high density regions.  
#3.) Save Leaf Clumps.  (Probably with Britton's tool. )
# .) For each clump
#4.)   Find clump indices from peaks
#5.)   Find data on prior snapshots
#6.)     -- Projections for simple vis
#7.)     -- Fake herschel, alma maps; local (region) imges, full cloud images.
#8.)     -- Physical Properties.

def clump_finder(ds, loc = None, width = None):
    if loc is None:
        value, loc = ds.find_max('density')
    if width is None:
        width = (0.05,'code_length')
    sphere = ds.sphere(loc,width)
    master_clump = Clump(sphere,"density")
    master_clump.add_validator("min_cells", 20)
    #master_clump.add_validator("gravitationally_bound", use_particles=False, use_thermal_energy=False)
    c_min = sphere["gas", "density"].min()
    c_max = sphere["gas", "density"].max()
    step = 100
    find_clumps(master_clump, c_min, c_max, step)
    return sphere, master_clump, loc

if 'u05' not in dir():
    u05 = taxi.taxi('u05')
    car = u05

#
# Get clumps.
#
if 0:
    car.fill(125)
    ds = car.ds
    sphere, master_clump, loc = clump_finder(ds,width=(0.05,'code_length'), loc = ds.arr([ 0.03613281,  0.79589844,  0.03027344], 'code_length'))
    proj = ds.proj('density',2,data_source=sphere,center=loc)
    leaf_clumps = get_lowest_clumps(master_clump) #if both min_cells and grav_bound are used, this is empty.
    pw = proj.to_pw(center = loc, width = (0.1,'code_length'))
    pw.annotate_clumps(clump_list_sort(leaf_clumps))
    pw.save('clump_testx')
elif 0:
    car.fill(125)
    ds = car.ds
    value, loc = ds.find_max('density')
    min_dx = ds.index.get_smallest_dx()
    region = ds.region(center=loc,left_edge = loc-1.5*min_dx, right_edge = loc+1.5*min_dx)
    indices = region['particle_index'].astype('int64')
    leaf_indices=[indices]
elif 1:
    car.fill(125)
    ds = car.ds
    keepers = [0,1,8,10,11,12,67,64,61, 201, 125, 306]
    colors  = ['r','g','b','c','m','y','r','g','b','c','m','y']
    peak_list = fPickle.load('u05_0125_peaklist.pickle')
    #indices_late,xpos_late,ypos_late,zpos_late  = {},{},{},{}
    leaf_indices={}
    """pull out indices for keepers.  Store in a dict"""
    for pind in keepers:
        peak = peak_list[pind]
        max_dx = ds.domain_width/ds.domain_dimensions
        min_dx = ds.index.get_smallest_dx().to_ndarray()
        locnd = peak # peak.to_ndarray()
        region = ds.region(center=peak,left_edge = peak-1.5*min_dx, right_edge = peak+1.5*min_dx)
        #indices_late[pind],xpos_late[pind],ypos_late[pind],zpos_late[pind] = clump_particles.particles_from_clump(region)
        leaf_indices[pind] = region['particle_index'].astype('int64')

#
# pre-images
#
frames_to_plot = [0,10,20 ,30,40,50,60,70,80,90,100,110,120,125][::2] #[10, 20, 30]:
cores_to_plot = keepers
cores_to_plot = [10,201]
if 1:
    for frame in frames_to_plot:
        car.fill(frame)
        ds  = car.ds
        proj_full = ds.proj('density',0,center='c')
        pw_full = proj_full.to_pw(center = 'c',width=(1.0,'code_length'))
        pw_full.set_cmap('density','gray')

        for nc,ic in enumerate(cores_to_plot):
            indices = leaf_indices[ic]
            pw_full.annotate_select_particles(1.0, col=colors[nc], indices=indices)
        print pw_full.save('u05_peak_thing_tmp3_%04d'%frame)

#
# early values, short version
#

if 0:
    for frame in [0]:
        car.fill(frame)
        ds  = car.ds
        for nc,pind in enumerate(cores_to_plot):
            ad = ds.all_data()
            indices = leaf_indices[pind]
            mask_to_get = np.zeros(indices.shape, dtype='int32') #this is necessary.  Factor of 2 faster to store mask_to_get
            t0=time.time()
            ad.set_field_parameter('indices_late',indices)
            ad.set_field_parameter('mask_to_get',mask_to_get)
            ad.set_field_parameter('timer',[t0])
            deposit_field = ad[("deposit",'deposit_target_particles')] #I believe this is necessary to build the field initially
            cut_region  = ad.cut_region(['obj["deposit","deposit_target_particles"] > 0'])
            mass = ( cut_region['density']*cut_region['cell_volume']).sum()
            print mass
            del ad

#
# long version.
#

if 0:
    if 'tdep' not in dir() or True:
        tdep = {}
    if 'qdict' not in dir() or True:
        qdict = {}
        cdict = {}
    #fields_to_plot = ['avg_density', 'kinetic_energy','binding_energy', 'avg_vx','avg_vy','avg_vz','mtotal','magnetic_energy']
    fields_to_plot = ['avg_density']
    times = []
    cycles=[]
    for field in fields_to_plot:
        if not qdict.has_key(field): qdict[field] = {}
        for c in cores_to_plot:
            qdict[field][c] = []
    for c in cores_to_plot:
        cdict[c]=clump_properties.series()
    n_cores=0
    for frame in frames_to_plot:
        #scratchdir = '/scratch1/dcollins/Paper19/u05-r4-l4-128'
        #scratchdir = '/work/00369/tg456484/maverick/Paper19/u05-r4-l4-128'
        #fname = '%s/DD%04d/data%04d'%(scratchdir,frame,frame)
        #ds = yt.load(fname)
        car.fill(frame)
        ds=car.ds
        times.append(ds['InitialTime'])
        cycles.append(ds['InitialCycleNumber'])
        if 'cut_region' in dir():
            del cut_region
        for n, pind in enumerate(cores_to_plot):
            if 'ad' in dir():
                del ad
            ad = ds.all_data()
            n_cores+=1
            mask_to_get = np.zeros(leaf_indices[pind].shape, dtype='int32') #this is necessary.  Factor of 2 faster to store mask_to_get
            t0=time.time()
            ad.set_field_parameter('indices_late',leaf_indices[pind])
            ad.set_field_parameter('mask_to_get',mask_to_get)
            ad.set_field_parameter('timer',[t0])
            #The dictionary behavior has not been tested
            tdep[pind] = ad[("deposit",'deposit_target_particles')]
            cut_region  = ad.cut_region(['obj["deposit","deposit_target_particles"] > 0'])
            #qdict['avg_density'][pind].append( (cut_region['density'].size))
            #print frame, pind
            #continue
            cdict[pind].append(clump_properties.clump_stuff(cut_region,ds,clump_2d_id = pind))
            if 0:
                print "energies on ", pind
                mtotal = cut_region['cell_mass'].sum()
                qdict['mtotal'][pind].append( mtotal)
                cell_mass = cut_region['cell_mass']
                cell_volume = cut_region['cell_volume']
                vx = cut_region['x-velocity']; vy = cut_region['y-velocity']; vz = cut_region['z-velocity']
                qdict['magnetic_energy'][pind].append( (cell_volume*(cut_region['Bx']**2+cut_region['By']**2 + cut_region['Bz']**2)).sum())
                vxbar=(vx*cell_mass).sum()/mtotal
                vybar=(vy*cell_mass).sum()/mtotal
                vzbar=(vz*cell_mass).sum()/mtotal
                qdict['avg_vx'][pind].append(vxbar)
                qdict['avg_vy'][pind].append(vybar)
                qdict['avg_vz'][pind].append(vzbar)
                qdict['avg_density'][pind].append( mtotal/cell_volume.sum())
                qdict['kinetic_energy'][pind].append( ( 0.5*cell_mass*( (vx-vxbar)**2 + (vy-vybar)**2 + (vz-vzbar)**2 )).sum())
                truncate = False
                if 0:
                  qdict['binding_energy'][pind] =ds['GravitationalConstant']* FindBindingEnergy(cut_region["gas", "cell_mass"].in_units('code_mass'),
                                                            cut_region["index", "x"].in_units('code_length'),
                                                            cut_region["index", "y"].in_units('code_length'),
                                                            cut_region["index", "z"].in_units('code_length'),
                                                            truncate, 1.0)
            del tdep[pind]

if 1:
    simname = 'test6'
    rmap = rainbow_map(n_cores)
    fields_to_plot = ['ge/ke','Alpha', 'Mass','Kinetic','Gravitational','Magnetic']
    fields_to_name = ['ge_ke','Alpha', 'Mass','ke','be','ge']
    cores_to_plot = [keepers[3], keepers[9]]
    for nf, field in enumerate(fields_to_plot):
        plt.clf()
        for nc, core in enumerate(cores_to_plot):
            if len(cdict[core][field]) == len(times):
                color = rmap(nc)
                color = colors[nc]
                plt.plot( times, cdict[core][field], color=color)
                plt.xlabel('t')
                plt.ylabel(field)
            else:
                print "Warning: error with core",core,"and field",field
        outname = 'p19_coresb_%s_%s.pdf'%(simname,fields_to_name[nf])
        plt.savefig(outname)
        print outname


