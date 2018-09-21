

frame= 1
d1 = '/scratch/00369/tg456484/Paper42_NewAK/dq04_m9_grav_512_L4_J4/DD%04d'%frame
d2 = '/scratch/00369/tg456484/Paper42_NewAK/dq04_m9_grav_512_L4_J4/RootGrids'

dsname1 = '%s/data%04d'%(d1,frame)
ds = yt.load(dsname1)
#fname1 = '%s/data%04d.cpu0000'%(d1,frame)
fname1 = ds.index.grids[0].filename

#RootOnly_broken_512_0001
fname2 = '%s/RootOnly_broken_512_%04d'%(d2,frame)
grid_id = 1
fptr_direct = h5py.File(fname1,'r')
fptr_cube   = h5py.File(fname2,'r')
field = 'z-velocity'
s1 = fptr_direct['Grid%08d'%grid_id][field][:]
s2 = fptr_cube[field][:]
fptr_direct.close()
fptr_cube.close()
print np.abs(s2[:64,:64,:64] - s1).sum()
print np.abs(s2[:64,:64,:64] + s1).sum()
