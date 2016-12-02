

setname = '/scratch1/dcollins/Paper42_new_turb/aq20_64_SRGIO_check/RootGrids/Root_x-velocity_y-velocity_z-velocity_64_test_0000'

fptr = h5py.File(setname)

try:
    vx = fptr['x-velocity'][:]
    vy = fptr['y-velocity'][:]
    vz = fptr['z-velocity'][:]

    sl1 = slice(
except:
    raise
finally:
    fptr.close()

