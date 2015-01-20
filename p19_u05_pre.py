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
scratchdir = '/work/00369/tg456484/maverick/Paper19/u05-r4-l4-128'
fname = '%s/DD%04d/data%04d'%(scratchdir,frame,frame)

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
        scratchdir = '/work/00369/tg456484/maverick/Paper19/u05-r4-l4-128'
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

    for field in ['avg_density', 'kinetic_energy','binding_energy', 'avg_vx','avg_vy','avg_vz','mtotal','magnetic_energy']:
        if not qdict.has_key(field): qdict[field] = {}

    for frame in [10]:
        scratchdir = '/scratch1/dcollins/Paper19/u05-r4-l4-128'
        scratchdir = '/work/00369/tg456484/maverick/Paper19/u05-r4-l4-128'
        fname = '%s/DD%04d/data%04d'%(scratchdir,frame,frame)
        ds = yt.load(fname)
        for n, pind in enumerate(keepers):
            ad = ds.all_data()
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
            qdict['mtotal'][pind] = mtotal
            cell_mass = cut_region['cell_mass']
            cell_volume = cut_region['cell_volume']
            vx = cut_region['x-velocity']; vy = cut_region['y-velocity']; vz = cut_region['z-velocity']
            qdict['magnetic_energy'][pind] = 0.5*(cell_volume*(cut_region['Bx']**2+cut_region['By']**2 + cut_region['Bz']**2)).sum()
            qdict['avg_vx'][pind] = (vx*cell_mass).sum()/mtotal
            qdict['avg_vy'][pind] = (vy*cell_mass).sum()/mtotal
            qdict['avg_vz'][pind] = (vz*cell_mass).sum()/mtotal
            qdict['avg_density'][pind] = mtotal/cell_volume.sum()
            qdict['kinetic_energy'][pind] =( 0.5*cell_mass*( (vx-qdict['avg_vx'][pind])**2 + (vy-qdict['avg_vy'][pind])**2 + (vz-qdict['avg_vz'][pind])**2 )).sum()
            truncate = False
            qdict['binding_energy'][pind] = FindBindingEnergy(cut_region["gas", "cell_mass"].in_units('code_mass'),
                                                        cut_region["index", "x"].in_units('code_length'),
                                                        cut_region["index", "y"].in_units('code_length'),
                                                        cut_region["index", "z"].in_units('code_length'),
                                                        truncate, 1.0)

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

