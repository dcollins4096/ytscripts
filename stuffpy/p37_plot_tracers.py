import clump_particles
d1 = '/Users/dcollins/scratch/Paper36_TracerTests/AddPost/b04_slab_fast'
frame = 49
setname = '%s/DD%04d/data%04d'%(d1,frame,frame)
ds=yt.load(setname)
center = [0.5]*3
sph = ds.sphere(center,0.1)
indices=sphere['particle_index']
proj = ds.proj('density',0,data_source=sph,center=center)
pw = proj.to_pw(center=center)
pw.set_cmap('density','gray')
pw.annotate_dave_particles(1,col='r',indices=indices)
print pw.save('b04_sph_%04d'%frame)
