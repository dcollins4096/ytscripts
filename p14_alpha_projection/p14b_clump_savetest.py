

from yt.analysis_modules.ppv_cube.api import PPVCube
import yt.units as u
import cPickle

if 1:
    nx,ny,nz = (256,256,256) # domain dimensions
    R = 10. # outer radius of disk, kpc
    r_0 = 3. # scale radius, kpc
    beta = 1.4 # for the tangential velocity profile
    alpha = -1. # for the radial density profile
    x, y = np.mgrid[-R:R:nx*1j,-R:R:ny*1j] # cartesian coordinates of x-y plane of disk
    r = np.sqrt(x*x+y*y) # polar coordinates
    theta = np.arctan2(y, x) # polar coordinates
    dens = np.zeros((nx,ny,nz))
    dens[:,:,nz/2-3:nz/2+3] = (r**alpha).reshape(nx,ny,1) # the density profile of the disk
    temp = np.zeros((nx,ny,nz))
    temp[:,:,nz/2-3:nz/2+3] = 1.0e5 # Isothermal
    vel_theta = 100.*r/(1.+(r/r_0)**beta) # the azimuthal velocity profile of the disk
    velx = np.zeros((nx,ny,nz))
    vely = np.zeros((nx,ny,nz))
    velx[:,:,nz/2-3:nz/2+3] = (-vel_theta*np.sin(theta)).reshape(nx,ny,1) # convert polar to cartesian
    vely[:,:,nz/2-3:nz/2+3] = (vel_theta*np.cos(theta)).reshape(nx,ny,1) # convert polar to cartesian
    dens[r > R] = 0.0
    temp[r > R] = 0.0
    velx[r > R] = 0.0
    vely[r > R] = 0.0
    data = {}
    data["density"] = (dens,"g/cm**3")
    data["temperature"] = (temp, "K")
    data["velocity_x"] = (velx, "km/s")
    data["velocity_y"] = (vely, "km/s")
    data["velocity_z"] = (np.zeros((nx,ny,nz)), "km/s") # zero velocity in the z-direction
    bbox = np.array([[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]]) # bbox of width 1 on a side with center (0,0,0)
    ds_round = yt.load_uniform_grid(data, (nx,ny,nz), length_unit=(2*R,"kpc"), nprocs=1, bbox=bbox)

    master_clump = Clump(ds_round.all_data(),'density')
    find_clumps(master_clump,16,300,2)
    leaf_clumps = get_lowest_clumps(master_clump) #if both min_cells and grav_bound are used, this is empty.

if 0:
    pickle_name = 'test_sphere_clump.pickle'
    fptr = open(pickle_name, 'wb')
    #Works and writes a set
    cPickle.dump(leaf_clumps[0].data,fptr, protocol=-1)
    fptr.close()


if 1:
    leaf_clumps[0].data.write_out('test_sphere_clumps.h5', fields=['density','velocity_x'])


if 0:
    pickle_name = 'test_sphere_clump.pickle'
    fptr = open(pickle_name, 'rb')
    new_clump=cPickle.load(fptr)
    fptr.close()
    #THis works if you run in the same python interpreter instance,
    #but fails if you launch a new python.
    print new_clump[1]['density']

