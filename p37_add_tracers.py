import re

#proj = yt.ProjectionPlot(ds,0,'density')
##proj.annotate_grids()
#proj.set_cmap('density','gray')
#proj.annotate_particles(1,col='r')
#print proj.save('derp')
#


class fake_grid_for_pointers():    
    def __init__(self):
        self.first_stuff=[]
        self.NumberOfParticles = None
        self.ParticleFileName = None
        self.BaryonFileName = None
        self.GravityBoundaryType = ''
        self.Pointers = []

def parse_hierarchy_for_particles(ds_name):
    #ds_name = get_ds_name(directory,frame)
    hname = "%s.hierarchy"%ds_name
    hptr = open(hname)
    grid_re = re.compile(r'^Grid = ([\d].*)')
    fake_grid_list = []
    try:
        for line in hptr:
            match = grid_re.match(line)
            if match:
                gnum= int(match.groups()[0])
                fake_grid_list.append(fake_grid_for_pointers())
                g_index = gnum-1
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("Task"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("GridRank"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("GridDimension"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("GridStartIndex"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("GridEndIndex"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("GridLeftEdge"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("GridRightEdg"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("Time"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("SubgridsAreStatic"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("NumberOfBaryonFields"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("FieldType"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("BaryonFileName"):
                fake_grid_list[g_index].first_stuff.append(line)
                fake_grid_list[g_index].BaryonFileName = line.split("=")[1]
            elif line.startswith("CourantSafetyNumber"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("PPMFlatteningParameter"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("PPMDiffusionParameter"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("PPMSteepeningParameter"):
                fake_grid_list[g_index].first_stuff.append(line)
            elif line.startswith("NumberOfParticles"):
                fake_grid_list[g_index].NumberOfPartiels = int(line.split("=")[1])
            elif line.startswith("ParticleFileName"):
                fake_grid_list[g_index].ParticleFileName = line.split("=")[1].strip()
            elif line.startswith("GravityBoundaryType"):
                fake_grid_list[g_index].GravityBoundaryType = line
            elif line.startswith("Pointer"):
                fake_grid_list[g_index].Pointers.append(line)
            elif len( line.strip() ): #it's not whitespace
                raise Exception("I don't recognize this line, thats super bad. \t\n%s\n"%line)

    except:
        raise
    finally:
        hptr.close()
    return fake_grid_list

def add_particles(ds, setname , outdir,outnumber=0, method=0):
    """adds a bunch of particles to ds.
    method = 0: one particle at each zone center."""
    fake_grid_list = parse_hierarchy_for_particles(setname)
    out_basename = outdir + "/DD%04d"%(outnumber)
    if glob.glob(out_basename) == []:
        os.mkdir(out_basename)

    out_hierarchy_name = outdir + "/DD%04d/%s%04d.hierarchy"%(outnumber,'data', outnumber)
    out_hierarchy_fptr = open(out_hierarchy_name,'w')
    particle_count_list = np.zeros(ds.index.grids.size)
    tracer_mass = 1e-10
    last_index = 0
    particle_fields = ['particle_index', 'particle_mass', 'particle_type',\
                       'particle_position_x', 'particle_position_y', 'particle_position_z', \
                       'particle_velocity_x', 'particle_velocity_y', 'particle_velocity_z']
    for grid_index, g in enumerate(ds.index.grids): 
        print "adding to grid",g,  'of', len(ds.index.grids)
        """identify the particles, make the lists"""
        list_of_lists = [[] for p in particle_fields]
        particle_dict = dict(zip(particle_fields,list_of_lists))
        density = g['density'].in_units('code_mass/code_length**3')#.flatten()
        x = g['x'].in_units('code_length')#flatten()
        y = g['y'].in_units('code_length')#flatten()
        z = g['z'].in_units('code_length')#flatten()
        vgx = g['x-velocity'].in_units('code_length/code_time')#flatten()
        vgy = g['y-velocity'].in_units('code_length/code_time')#flatten()
        vgz = g['z-velocity'].in_units('code_length/code_time')#flatten()
        
        this_number = particle_count_list[grid_index] = g.child_mask.sum()
        particle_dict['particle_type'] = np.zeros(this_number)+3
        particle_dict['particle_mass'] = np.zeros(this_number)+tracer_mass
        particle_dict['particle_index'] = np.arange(last_index, last_index+this_number)
        last_index += this_number
        particle_dict['particle_position_x'] = x[g.child_mask].flatten()
        particle_dict['particle_position_y'] = y[g.child_mask].flatten()
        particle_dict['particle_position_z'] = z[g.child_mask].flatten()
        particle_dict['particle_velocity_x'] = vgx[g.child_mask].flatten()
        particle_dict['particle_velocity_y'] = vgy[g.child_mask].flatten()
        particle_dict['particle_velocity_z'] = vgz[g.child_mask].flatten()
        #mask = g.child_index_mask.flatten()
#       got_some = 0
#       for n in range(density.size):
#           if mask[n] < 0: # and got_some < 3:
#               got_some += 1
#               particle_dict['particle_type'].append(3)
#               particle_dict['particle_mass'].append(tracer_mass)
#               particle_dict['particle_index'].append(last_index )
#               last_index += 1
#               particle_dict['particle_position_x'].append(x[n])
#               particle_dict['particle_position_y'].append(y[n])
#               particle_dict['particle_position_z'].append(z[n])
#               particle_dict['particle_velocity_x'].append(vgx[n])
#               particle_dict['particle_velocity_y'].append(vgy[n])
#               particle_dict['particle_velocity_z'].append(vgz[n])
        del density , x, y, z, vgx, vgy, vgz
#        particle_count_list[grid_index] = len(particle_dict['particle_type'])
        if fake_grid_list[grid_index].NumberOfParticles is not None:
            particle_count_list[grid_index]+= fake_grid_list[grid_index].NumberOfParticles
        
        """Write the data files"""
        out_file_name = "%s/%s"%(out_basename , g.filename.split("/")[-1])
        print "writing", out_file_name
        in_cpu = h5py.File(g.filename,'r')
        in_group = in_cpu['Grid%08d'%g.id]
        out_cpu = h5py.File(out_file_name,'a')
        out_group = out_cpu.require_group( 'Grid%08d'%g.id )
        for in_field in in_group:
            if in_field in particle_fields:
                """Concatenate the existing particles"""
                particle_dict[in_field].append(in_group[in_field])
            else:
                out_group.require_dataset(name=in_field,shape=in_group[in_field].shape,dtype=in_group[in_field].dtype)
                out_group[in_field][:] = in_group[in_field][:]
            
        for field in particle_fields:
            dtype = np.dtype('f8')
            if field in ['particle_index','particle_type']:
                dtype = np.dtype('i8')
            out_group.require_dataset(name=field,shape=(int(particle_count_list[grid_index]),) ,dtype=dtype)
            out_group[field][:] = particle_dict[field]

        out_cpu.close()
        in_cpu.close()
        for field in particle_fields:
            del particle_dict[field]
        
        

        """write the hierarchy files"""
        for line in fake_grid_list[grid_index].first_stuff:
            out_hierarchy_fptr.write(line)
        NumberOfParticles = particle_count_list[grid_index]
        out_hierarchy_fptr.write("NumberOfParticles = %d\n"%(NumberOfParticles))
        if NumberOfParticles > 0 :
            out_hierarchy_fptr.write("ParticleFileName = %s"%fake_grid_list[grid_index].BaryonFileName) 
        out_hierarchy_fptr.write(fake_grid_list[grid_index].GravityBoundaryType)
        for line in fake_grid_list[grid_index].Pointers:
            out_hierarchy_fptr.write(line)
        out_hierarchy_fptr.write("\n")

        
        
    out_hierarchy_fptr.close()
    """write the parameter files"""
    out_ds_name = outdir + "/DD%04d/%s%04d"%(outnumber,'data',outnumber)
    out_ds = open(out_ds_name,'w')
    in_ds = open(setname,'r')
    for line in in_ds:
        if line.startswith('TracerParticleOn'):
            out_ds.write('TracerParticleOn = 1\n')
        elif line.startswith('NumberOfParticles'):
            out_ds.write('NumberOfParticles      = %d (do not modify)\n'%particle_count_list.sum()) #haha, do not modify.
        elif line.startswith('DataDumpNumber'):
            out_ds.write('DataDumpNumber = %d\n'%(outnumber+1))
        elif line.startswith('StopCycle'):
            StopCycle = int(line.split('=')[1])
            out_ds.write('StopCycle = %d\n'%(StopCycle+3))
        elif line.startswith('CycleSkipDataDump'):
            CycleSkipDataDump = int(line.split('=')[1])
            out_ds.write('CycleSkipDataDump = %d\n'%(1))

        else:
            out_ds.write(line)
    out_ds.close()
    in_ds.close()
    """get all the damned anscilary files.  Should be done with copies."""
    extra_files=['.boundary','.configure','.boundary.hdf','.memorymap']
    for fl in extra_files:
        source_file = "%s%s"%(setname,fl)
        dest_file = "%s%s"%(out_ds_name,fl)
        shutil.copy(source_file,dest_file)

def copy_driving(ds1, ds2, fields=["%s-acceleration"%s for s in "xyz"]):
    """adds the [xyz]-acceleration from ds1 to ds2"""

    for grid_index, g1 in enumerate(ds1.index.grids): 
        g2 = ds2.index.grids[grid_index]
        if np.abs(g1.LeftEdge - g2.LeftEdge).sum() > g2.dds.min():
            raise CopyError("Left Edge")
        if np.abs(g1.RightEdge - g2.RightEdge).sum() > g2.dds.min():
            raise CopyError("Left Edge")
        if np.abs( g1.dds.in_units('code_length') - g2.dds.in_units('code_length')).sum() > g2.dds.in_units('code_length').min()*1e-7:
            raise CopyError("Cell Width")
        if np.abs( g1.ActiveDimensions  - g2.ActiveDimensions ).sum() > 0:
            raise CopyError("Nzones")

        in_cpu = h5py.File(g1.filename,'r')
        in_group = in_cpu['Grid%08d'%g1.id]
        out_cpu = h5py.File(g2.filename,'r+')
        out_group = out_cpu.require_group( 'Grid%08d'%g2.id )
        for field in fields:
            out_group.require_dataset(name=field,shape=in_group[field].shape, dtype=in_group[field].dtype)
            out_group[field][:] = in_group[field][:]
        out_cpu.close()
        in_cpu.close()


    setname = ds2.fullpath+"/"+str(ds2)
    fake_grid_list = parse_hierarchy_for_particles(setname)
    out_hierarchy_name = setname + ".hierarchy.new"
    in_hierarchy_name = setname + ".hierarchy"
    out_hierarchy_fptr = open(out_hierarchy_name,'w')
    in_hierarchy_fptr = open(in_hierarchy_name,'r')
    for line in in_hierarchy_fptr:
        if line.startswith("FieldType"):
            out_hierarchy_fptr.write( line[:-1] + " 55 56 57\n")
        elif line.startswith("NumberOfBaryonFields"):
            N = int(line.split(" ")[-1]) +3
            out_hierarchy_fptr.write("NumberOfBaryonFields = %d\n"%N)


        else:
            out_hierarchy_fptr.write(line)
    out_hierarchy_fptr.close()







ds1 = yt.load('/scratch1/dcollins/Paper48_DrivingFiller/Breeder8/DD0041/data0041')
ds2 = yt.load('/scratch1/dcollins/Paper48_DrivingFiller/DD0000/data0000')
copy_driving(ds1,ds2)


#dirname = '/scratch/00369/tg456484/Paper37_Restart/a00_ics'
#outdir  = '/scratch/00369/tg456484/Paper37_Restart/a01_with_uniform_particles'
#dirname = '/Users/dcollins/scratch/Paper36_TracerTests/AddPost/b01_slab'
#outdir = '/Users/dcollins/scratch/Paper36_TracerTests/AddPost/b03_slab_add_particles'
#outdir = '/Users/dcollins/scratch/Paper36_TracerTests/AddPost/b04_slab_fast'
#frame = 60;setname = '%s/RS%04d/restart%04d'%(dirname,frame,frame) 

if 0:
    """works."""
    dirname = '/scratch1/dcollins/Paper36_tracertests/AddPost/c02_small_sphere_noparticles'
    outdir = '/scratch1/dcollins/Paper36_tracertests/AddPost/c03_c02_added'
    frame = 0;setname = '%s/DD%04d/data%04d'%(dirname,frame,frame)
    ds = yt.load(setname)
    add_particles(ds,setname,outdir)
if 0:
    """works."""
    dirname = '/scratch/00369/tg456484/Paper42_NewAK/bq42_m9_drive0.5_512_p59_ppm/'
    outdir = '/scratch/00369/tg456484/Paper42_NewAK/bq42_m9_drive0.5_512_p59_ppm/'
    frame = 40;setname = '%s/IC%04d/data%04d'%(dirname,frame,frame)
    ds = yt.load(setname)
    add_particles(ds,setname,outdir)
