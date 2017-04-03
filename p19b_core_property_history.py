if 'ef' not in dir():
    execfile('go')
from yt.analysis_modules.level_sets.api import * #for clumps
import pyximport; pyximport.install()
import particle_ops
import random
import clump_particles
from yt.utilities.data_point_utilities import FindBindingEnergy
from yt.utilities.physical_constants import \
            gravitational_constant_cgs as G
reload(clump_particles)
ef('particle_selector.py')

frame = 125
scratchdir = '/scratch1/dcollins/Paper19/u05-r4-l4-128'
#scratchdir = '/work/00369/tg456484/maverick/Paper19/u05-r4-l4-128'
fname = '%s/DD%04d/data%04d'%(scratchdir,frame,frame)
simname = 'u05'

if 'ds' not in dir():
    ds = yt.load(fname)

if 'peak_list' not in dir():
    peak_list = fPickle.load('u05_0125_peaklist.pickle')

keepers = [0,1,8,10,11,12,67,64,61, 201, 125, 306]
colors  = ['r','g','b','c','m','y','r','g','b','c','m','y']
if 'indices_late' not in dir():
    indices_late,xpos_late,ypos_late,zpos_late  = {},{},{},{}
    """pull out indices for keepers.  Store in a dict"""
    for pind in keepers:
        peak = peak_list[pind]
        max_dx = ds.domain_width/ds.domain_dimensions
        min_dx = ds.index.get_smallest_dx().to_ndarray()
        locnd = peak # peak.to_ndarray()
        region = ds.region(center=peak,left_edge = peak-1.5*min_dx, right_edge = peak+1.5*min_dx)
        #indices_late[pind],xpos_late[pind],ypos_late[pind],zpos_late[pind] = clump_particles.particles_from_clump(region)
        indices_late[pind] = region['particle_index'].astype('int64')

if 0:
    """plot the particles for a sequence of frames."""
    ax = 1
    for frame in [0,10,20,30,40,50,60,70,80,90,100,110,120,125]:
        scratchdir = '/scratch1/dcollins/Paper19/u05-r4-l4-128'
        #scratchdir = '/work/00369/tg456484/maverick/Paper19/u05-r4-l4-128'
        fname = '%s/DD%04d/data%04d'%(scratchdir,frame,frame)
        ds = yt.load(fname)
        proj_full = ds.proj('density',ax, center = 'c' ) #width = (1.0,'code_length'))
        pw_full = proj_full.to_pw(center = 'c',width=(1.0,'code_length'))
        pw_full.set_cmap('density','gray')
        for n, pind in enumerate(keepers):
            pw_full.annotate_dave_particles(1.0, col=colors[n], indices=indices_late[pind])
            if frame == 125:
                pw_full.annotate_point(peak_list[pind],pind,text_args={'color':colors[n]})
        print pw_full.save('u05_%04d_keepers1'%frame)

if 1:
    if 'tdep' not in dir():
        tdep = {}
    if 'qdict' not in dir():
        qdict = {}
    def flatten(dic,field):
        fld = dic[field]
        keys = sorted(fld.keys())
        out = []
        for k in keys:
            try:
                out.append( fld[k].item() )
            except:
                out.append( fld[k] )

        return nar(out)

    fields_to_plot = ['avg_density', 'kinetic_energy','binding_energy', 'avg_vx','avg_vy','avg_vz','mtotal','magnetic_energy']
    times = []
    cycles=[]
    cores_to_plot = keepers #[0:3]
    for field in fields_to_plot:
        if not qdict.has_key(field): qdict[field] = {}
        for c in cores_to_plot:
            qdict[field][c] = []

    n_cores=0
    for frame in [0,10,20,30,40,50,60,70,80,90,100,110,120,125]:#[10, 20, 30]:
        scratchdir = '/scratch1/dcollins/Paper19/u05-r4-l4-128'
        #scratchdir = '/work/00369/tg456484/maverick/Paper19/u05-r4-l4-128'
        fname = '%s/DD%04d/data%04d'%(scratchdir,frame,frame)
        ds = yt.load(fname)
        ad = ds.all_data()
        times.append(ds['InitialTime'])
        cycles.append(ds['InitialCycleNumber'])
        for n, pind in enumerate(cores_to_plot):
            n_cores+=1
            mask_to_get = np.zeros(indices_late[pind].shape, dtype='int32') #this is necessary.  Factor of 2 faster to store mask_to_get
            t0=time.time()
            ad.set_field_parameter('indices_late',indices_late[pind])
            ad.set_field_parameter('mask_to_get',mask_to_get)
            ad.set_field_parameter('timer',[t0])
            #The dictionary behavior has not been tested.
            tdep[pind] = ad[("deposit",'deposit_target_particles')]
            cut_region  = ad.cut_region(['obj["deposit","deposit_target_particles"] > 0'])
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
    rmap = rainbow_map(n_cores)
    for nf, field in enumerate(fields_to_plot):
        plt.clf()
        for nc, core in enumerate(cores_to_plot):
            if len(qdict[field][core]) == len(times):
                plt.plot( times, qdict[field][core], color=rmap(nc))
                plt.xlabel('t')
                plt.ylabel(field)
            else:
                print "Warning: error with core",core,"and field",field
        outname = 'p19_cores_%s_%s.pdf'%(simname,field)
        plt.savefig(outname)
        print outname


#try to get the preimage
if 0:
    #This does not work.
    master_clump = Clump(ad,'dp1')
    #master_clump = Clump(ad,("deposit",'deposit_target_particles_1'))
    #master_clump = Clump(ad,("deposit",'all_cic'))
    c_min = 0.9
    c_max = indices_late[pind].size
    step = c_max*1.01
    print "CLOWN make the clumps"
    find_clumps(master_clump, c_min, c_max, step)

