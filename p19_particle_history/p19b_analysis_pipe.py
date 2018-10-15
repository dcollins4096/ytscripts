#execfile('go_lite')
execfile('go')
import taxi
from p19b_select_particle_callback import *
#import clump_particles #
import time
import clump_properties
reload(clump_properties)
import pyximport; pyximport.install()
import particle_grid_mask
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

def shift_particles(ds, position,shiftRight = True):
    """Shifts a periodically separated clump by the domain width.
    Looks for gaps in the positions larger than max('dx'), shifts one group
    to the right (left if shiftRight=Flase) to be spatially contiguous."""
    #max_dx may not be computed in the most efficient way.
    DomainLeft = ds.domain_left_edge
    DomainRight =  ds.domain_right_edge
    DomainWidth = DomainRight - DomainLeft
    shift = np.zeros(3)
    for i,axis in enumerate(['x','y','z']):
        if i != 1:
            continue

        dx = 'd'+axis
        nique = np.unique(position[i])
        nique.sort()
        print(nique)
        max_dx = ds.index.grids[0].dds.max().in_units('code_length')
        min_dx = ds.index.get_smallest_dx()

        #has to be close to the edges, or 'Periodic Wrap' isn't the problem.
        if np.abs(nique.max() - DomainRight[i]) > 3*max_dx:
            continue
        if np.abs(nique.min() - DomainLeft[i]) > 3*max_dx:
            continue
        delta_x = nique[1:] - nique[0:-1]
        max_delta_x = delta_x.max()
        break_index = np.where(delta_x == max_delta_x)
        print("SHIFT max %0.3f index %d"%(max_delta_x,break_index[0][0]))

        if max_delta_x > max_dx:
           break_x = nique[break_index[0]+1]
           print("SHIFT break x %0.3f"%break_x)
           if shiftRight:
              all_to_shift = np.where( position[i] <= break_x + min_dx )[0]
              position[i][all_to_shift] += DomainWidth[i]
              shift[i] = DomainWidth[i]
           else:
                all_to_shift = np.where( position[i] >= break_x - min_dx )[0]
                position[i][all_to_shift] -= DomainWidth[i]
                shift[i] = -DomainWidth[i]
                
    return all_to_shift

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

#if 'u05' not in dir():
#    u05 = taxi.taxi('u05')
#    car = u05
car = taxi.taxi('u05')
#
# Get clumps.
#
if 0:
    ds = car.load(157)
    sphere, master_clump, loc = clump_finder(ds,width=(0.05,'code_length'), loc = ds.arr([ 0.03613281,  0.79589844,  0.03027344], 'code_length'))
    proj = ds.proj('density',2,data_source=sphere,center=loc)
    leaf_clumps = get_lowest_clumps(master_clump) #if both min_cells and grav_bound are used, this is empty.
    pw = proj.to_pw(center = loc, width = (0.1,'code_length'))
    pw.annotate_clumps(clump_list_sort(leaf_clumps))
    pw.save('clump_testx')
elif 0:
    ds=car.load(125)
    value, loc = ds.find_max('density')
    min_dx = ds.index.get_smallest_dx()
    region = ds.region(center=loc,left_edge = loc-1.5*min_dx, right_edge = loc+1.5*min_dx)
    indices = region['particle_index'].astype('int64')
    leaf_indices=[indices]
elif 1:
    ds = car.load(125)
    keepers = [0,1,8,10,11,12,67,68,64,61, 201, 125, 306]
    colors  = ['r','g','b','c','m','y','r','g','b','c','m','y']
    peak_list = fPickle.load('u05_0125_peaklist.pickle')
    #indices_late,xpos_late,ypos_late,zpos_late  = {},{},{},{}
    leaf_indices={}
    """pull out indices for keepers.  Store in a dict"""
    for core_ind in keepers:
        peak = peak_list[core_ind]
        max_dx = ds.domain_width/ds.domain_dimensions
        min_dx = ds.index.get_smallest_dx().to_ndarray()
        locnd = peak # peak.to_ndarray()
        region = ds.region(center=peak,left_edge = peak-1.5*min_dx, right_edge = peak+1.5*min_dx)
        #indices_late[core_ind],xpos_late[core_ind],ypos_late[core_ind],zpos_late[core_ind] = clump_particles.particles_from_clump(region)
        leaf_indices[core_ind] = region['particle_index'].astype('int64')

#
# cherry pick for debugging.
#

frames_to_plot = [10] # range(0,155,20)+[155]# [0,10,20 ,30,40,50,60,70,80,90, 100 ,110,120,125] #[120] # [0,10,20 ,30,40,50,60,70,80,90,100,110,120,125][::2] #[10, 20, 30]:
#cores_to_plot = keepers
cores_to_plot = [8] # [0,1,8,68]

#
# pre-images
#
if 0:
    for frame in frames_to_plot:
        ds = car.load(frame)
        ad=ds.all_data()
        this_center=[0.5,0.0,0.5]
        proj_full = ds.proj('density',0,center=this_center)
        pw_full = proj_full.to_pw(center = this_center,width=(1.0,'code_length'),origin='domain')
        pw_full.set_cmap('density','gray')
        print pw_full.save('%s_peak_thing_core_center_scatter_%04d'%(car.name,frame))
        for nc,ic in enumerate(cores_to_plot):
            indices = leaf_indices[ic]
            pw_full.annotate_select_particles(1.0, col=colors[nc], indices=indices)
            pw_full.annotate_text(peak_list[ic],'%d'%ic,text_args={'color':'r'})
            these_pids = leaf_indices[ic].astype('int64').v
            #ef('tmp.py')
            my_indices = ad['particle_index'].astype('int64')
            mask_to_get_2 = np.zeros(indices.shape, dtype='int32')
            found_any, mask = particle_ops.mask_particles(
                these_pids, my_indices, mask_to_get_2)
            pos = ad['particle_position'][mask == 1]
            pos_1 = copy.copy(pos)
            xpos,ypos,zpos = pos[:,0], pos[:,1], pos[:,2]
            xpos_1,ypos_1,zpos_1 = pos_1[:,0], pos_1[:,1], pos_1[:,2]
            sha = shift_particles(ds,[xpos,ypos,zpos], shiftRight=False)
            xbar = xpos.mean()
            ybar = ypos.mean()
            zbar = zpos.mean()
            rmax = np.max( (xpos-xbar)**2+(ypos-ybar)**2+(zpos-zbar)**2)**0.5
            center_sph = [xbar,ybar,zbar]
            sph = ds.sphere(center_sph, rmax)
            pw_full.annotate_sphere(center_sph,rmax)
            pw_full.save('tmp.png')
            pw_full.plots['density'].axes.scatter(ypos,zpos)
            #pw_full.save('test5.png')
#           print("poo")
            #pw_full.plots['density'].axes.scatter(np.arange(0,1,0.1),np.arange(0,1,0.1))

        #pw_full.plots['density'].axes.scatter(np.arange(0,1,0.1),np.arange(0,1,0.1))
        print pw_full.save('%s_peak_thing_core_center_scatter_%04d'%(car.name,frame))
        plt.clf()
        plt.scatter(ypos,zpos)
        plt.savefig('test7.png')
        plt.clf()
        plt.scatter(ypos_1,zpos_1)
        plt.scatter(ypos_1[sha],zpos_1[sha],c='r')
        plt.savefig('test71.png')

#
# long version.
#
if 0:
    if 'tdep' not in dir() or True:
        tdep = {}
    if 'qdict' not in dir() or True:
        qdict = {}
        cdict = {}
        pdict = {}
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
        for n, core_ind in enumerate(cores_to_plot):
            if 'ad' in dir():
                del ad
            ad = ds.all_data()
            n_cores+=1
            if 0:
                """ Plotting of average clump properties."""
                mask_to_get = np.zeros(leaf_indices[core_ind].shape, dtype='int32') #this is necessary.  Factor of 2 faster to store mask_to_get
                t0=time.time()
                ad.set_field_parameter('indices_late',leaf_indices[core_ind])
                ad.set_field_parameter('mask_to_get',mask_to_get)
                ad.set_field_parameter('timer',[t0])
                #The dictionary behavior has not been tested
                tdep[core_ind] = ad[("deposit",'deposit_target_particles')]
                cut_region  = ad.cut_region(['obj["deposit","deposit_target_particles"] > 0'])
                print """ Look at the .ori verison.  Somethin brok."""
                #qdict['avg_density'][core_ind].append( (cut_region['density'].size))
                #print frame, core_ind
                #continue
                cdict[core_ind].append(clump_properties.clump_stuff(cut_region,ds,clump_2d_id = core_ind))
            if 0:
                """ a thing that didn't work."""
                particle_index = x + nx*(y + nx*z)
                cx,cy,cz = [(cut_region[ax]/min_dx).astype('int64') for ax in 'xyz']
                cx_index = cx + nx*(cy + nx*cz)
                c_argsort = np.argsort(cx_index)
                p_argsort = np.argsort(particle_index)
                psort = particle_index[p_argsort]
            if 1:
                #dan look at this bit here.
                these_pids = leaf_indices[core_ind].astype('int64').v
                #ef('tmp.py')
                my_indices = ad['particle_index'].astype('int64')
                mask_to_get_2 = np.zeros(these_pids.shape, dtype='int32')
                found_any, mask = particle_ops.mask_particles(
                    these_pids, my_indices, mask_to_get_2)
                pos = ad['particle_position'][mask == 1]
                xpos,ypos,zpos = pos[:,0], pos[:,1], pos[:,2]
                for grid in ds.index.grids[-1::-1]:
                    grid_selector = np.zeros([3,these_pids.size], dtype='int32')
                    particle_selector = np.zeros(these_pids.size, dtype='int32')
                    mask_to_get_3 = np.zeros(these_pids.shape, dtype='int32')
                    found_any_g, mask_g = particle_ops.mask_particles(
                        these_pids, grid['particle_index'].astype('int64'), mask_to_get_3)
                    if found_any_g:
                        particle_grid_mask.particle_grid_mask_go_i2( #this may be broken.
                            xpos,ypos,zpos, grid.LeftEdge, grid.dds, grid.ActiveDimensions,grid.child_mask, grid_selector,particle_selector)
                        particle_selector = particle_selector == 1
                        particles_in_grid = these_pids[particle_selector]
                        for p in particles_in_grid:
                            pdict[p] = pdict.get(p,{})
                            for k in ['times','cycles']:
                                if not pdict[p].has_key(k): pdict[p][k]=[]
                            pdict[p]['times'].append(times[-1])
                            pdict[p]['cycles'].append(cycles[-1])
                        for field in ['density']:
                            values = grid[field][[grid_selector[i][particle_selector] for i in [0,1,2]]]
                            for n,p in enumerate(particles_in_grid):
                                pdict[p][field] = pdict[p].get(field,[])
                                pdict[p][field].append(values[n])
            if 0:
                """straight and direct grabbing of values.  Not great."""
                print "energies on ", core_ind
                mtotal = cut_region['cell_mass'].sum()
                qdict['mtotal'][core_ind].append( mtotal)
                cell_mass = cut_region['cell_mass']
                cell_volume = cut_region['cell_volume']
                vx = cut_region['x-velocity']; vy = cut_region['y-velocity']; vz = cut_region['z-velocity']
                qdict['magnetic_energy'][core_ind].append( (cell_volume*(cut_region['Bx']**2+cut_region['By']**2 + cut_region['Bz']**2)).sum())
                vxbar=(vx*cell_mass).sum()/mtotal
                vybar=(vy*cell_mass).sum()/mtotal
                vzbar=(vz*cell_mass).sum()/mtotal
                qdict['avg_vx'][core_ind].append(vxbar)
                qdict['avg_vy'][core_ind].append(vybar)
                qdict['avg_vz'][core_ind].append(vzbar)
                qdict['avg_density'][core_ind].append( mtotal/cell_volume.sum())
                qdict['kinetic_energy'][core_ind].append( ( 0.5*cell_mass*( (vx-vxbar)**2 + (vy-vybar)**2 + (vz-vzbar)**2 )).sum())
                truncate = False
                if 0:
                  qdict['binding_energy'][core_ind] =ds['GravitationalConstant']* FindBindingEnergy(cut_region["gas", "cell_mass"].in_units('code_mass'),
                                                            cut_region["index", "x"].in_units('code_length'),
                                                            cut_region["index", "y"].in_units('code_length'),
                                                            cut_region["index", "z"].in_units('code_length'),
                                                            truncate, 1.0)
            if False and tdep.has_key(core_ind):
                del tdep[core_ind]
#ef('p19b_plot_pipe.py')
#
# early values, short version
#
if 0:
    for frame in [0]:
        car.fill(frame)
        ds  = car.ds
        for nc,core_ind in enumerate(cores_to_plot):
            ad = ds.all_data()
            indices = leaf_indices[core_ind]
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

