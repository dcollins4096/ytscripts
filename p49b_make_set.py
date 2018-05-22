if 'ef' not in dir():
    execfile('go')
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
d=1.
vx=vy=vz=0
Bx, By, Bz = 1.0,1.41421,0.5
P = 0.6
Gamma=1.6666666667
Egas=P/((Gamma-1)*d)
E = 0.5*(vx*vx+vy*vy+vz*vz)+0.5*(Bx*Bx+By*By+Bz*Bz)/d + Egas
right = p49_eigen.eigen(d, vx, vy, vz, Bx,By,Bz,E, Gamma=Gamma)
speeds = p49_eigen.speeds(d, vx, vy, vz, Bx,By,Bz,E, Gamma=Gamma)
EigenVector = np.zeros([7,7])
if 0:
    CLOWN_EigenVector = np.zeros([7,7])
    CLOWN_EigenVector[0, 0]=1.000000
    CLOWN_EigenVector[1, 0]=1.640000
    CLOWN_EigenVector[2, 0]=-1.000000
    CLOWN_EigenVector[3, 0]=0.707107
    CLOWN_EigenVector[4, 0]=0.707107
    CLOWN_EigenVector[5, 0]=0.707107
    CLOWN_EigenVector[6, 0]=0.707107
if 0:
    EigenVector[0,0]=0.894427
    EigenVector[1,0]=4.024912
    EigenVector[2,0]=-1.788851
    EigenVector[3,0]=0.843274
    EigenVector[4,0]=0.298143
    EigenVector[5,0]=1.686545
    EigenVector[6,0]=0.596285
if 0:
    EigenVector[0,0]= 0.892498 
    EigenVector[1,0]= 4.112684 
    EigenVector[2,0]= -1.777855 
    EigenVector[3,0]= 0.847106 
    EigenVector[4,0]= 0.299498 
    EigenVector[5,0]= 1.687435 
    EigenVector[6,0]= 0.596600 
if 1:
    EigenVector[0,4]= 0.898492 
    EigenVector[1,4]= 0.792409 
    EigenVector[2,4]= 0.441938 
    EigenVector[3,4]= 0.824454 
    EigenVector[4,4]= 0.291489 
    EigenVector[5,4]= -0.405521 
    EigenVector[6,4]= -0.143374 

for i,j in zip(right[:,4], EigenVector[:,4]):
    print(i,j, (i-j)/(i+j))

"""
shape = np.array([16,16,16])
dx = 1./shape
x,y,z=np.mgrid[0:1:dx[0],0:1:dx[1],0:1:dx[2]]
density=np.ones(shape)
vx = np.zeros(shape)
vy = np.zeros(shape)
vz = np.zeros(shape)

directory = '/Users/dcollins/scratch/Paper49b_play/r13d_left_fast_new'
if 0:
    r_circ = 0.2
    center=[0.5,0.5]
    r_vec = ((x-center[0])**2+(y-center[1])**2)**0.5
    mask = r_vec < r_circ
    vx[mask] = -0.1*(y[mask]-center[1])
    vy[mask] = 0.1*(x[mask]-center[0])
    density=1+x

    directory = '/Users/dcollins/scratch/Paper49b_play/r13b_density_test/'

    xx,yx,zx=np.mgrid[0:1+dx[0]:dx[0],0:1:dx[1],   0:1:dx[2]]
    xy,yy,zy=np.mgrid[0:1:dx[0],      0:1+dx[1]:dx[1],0:1:dx[2]]
    xz,yz,zz=np.mgrid[0:1:dx[0],      0:1:dx[1],   0:1+dx[2]:dx[2]]

    Bx = np.zeros_like(xx)
    Bx = yx
    By = np.zeros_like(xy)
    Bz = np.zeros_like(xz)

if 1:
    p49_eigen.eigen(1,0,0,0,0,0,0,1)


file_density = '%s/test_density'%directory
file_vx = '%s/test_vx'%directory
file_vy = '%s/test_vy'%directory
file_vz = '%s/test_vz'%directory
file_Bx = '%s/test_Bx'%directory
file_By = '%s/test_By'%directory
file_Bz = '%s/test_Bz'%directory

enzo_write.dump_h5(density,file_density)
enzo_write.dump_h5(vx,file_vx)
enzo_write.dump_h5(vy,file_vy)
enzo_write.dump_h5(vz,file_vz)

enzo_write.dump_h5(Bx,file_Bx)
enzo_write.dump_h5(By,file_By)
enzo_write.dump_h5(Bz,file_Bz)
"""
