ef('clump_particles.py')
basedir = '/scratch/00369/tg456484/Paper37_Restart/B02/512'
if 1:
    set50 = '%s/RS%04d/restart%04d'%(basedir,50,50)
    set53 = '%s/RS%04d/restart%04d'%(basedir,51,51)
    ds50 = yt.load(set50)
    ds53 = yt.load(set53)

frame=50
if 0:
    center = [0.5]*3
    width = 0.01
if 0:
    max_density, center = ds50.find_max('density')
if 1:
    #max_density, center = (34993852.0 g/cm**3, YTArray([ 0.59735107,  0.73248291,  0.96258545]) code_length)
    max_density, center = (34993852., np.array([ 0.59735107,  0.73248291,  0.96258545]) )
    width = 0.01

if 0:
    sph = ds50.sphere(center,width)
    indices50=sph['particle_index']
    proj = ds50.proj('density',0,data_source=sph,center=center)
    pw = proj.to_pw(center=center, width = 2*width)
    pw.set_cmap('density','gray')
    pw.annotate_dave_particles(1,col='r',indices=indices50)
    print pw.save('p37_r1_sph_%04d'%frame)

if 0:
    frame=53
    sph53 = ds53.sphere(center,width)
    indices53=sph53['particle_index']
    proj = ds53.proj('density',0,data_source=sph53,center=center)
    pw = proj.to_pw(center=center, width = 2*width)
    pw.set_cmap('density','gray')
    pw.annotate_dave_particles(1,col='r',indices=indices53)
    #pw.annotate_particles(1,col='r') #,indices=indices50)
    print pw.save('p37_r1b_sph_ind53_%04d'%frame)

if 1:
    """Are particles missing from upper levels?"""
    n_particles50 = np.zeros(6)
    n50 = len(ds50.index.grids)
    for g in ds50.index.grids:
        print g, "of" , n50
        n_particles50[g.Level] += g['particle_index'].size
        del g['particle_index']
    n_particles53 = np.zeros(6)
    n53 = len(ds53.index.grids)
    for g in ds53.index.grids:
        print g, "of" , n53
        n_particles53[g.Level] +=g['particle_index'].size
        #del g['particle_index']
    print n_particles50
    print n_particles53


