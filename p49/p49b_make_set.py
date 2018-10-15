if 'ef' not in dir():
    execfile('go')
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
class system():
    def __init__(self, d=1.0,vx=0.0,vy=0.0,vz=0.0,Bx=1.0,By=1.41421,Bz=0.5,Gamma=1.6666666667,E=None,P=None, write=True):
        #P = 0.6
        #Gamma=1.6666666667
        if P is not None:
            Egas=P/((Gamma-1)*d)
            E = 0.5*(vx*vx+vy*vy+vz*vz)+0.5*(Bx*Bx+By*By+Bz*Bz)/d + Egas
        self.quan={'d':d,'vx':vx,'vy':vy,'vz':vz,'Bx':Bx,'By':By,'Bz':Bz,'P':P,'Gamma':Gamma,'Egas':Egas,'E':E}
        for a,b in [['d','Density'],['vx','x-velocity'],['vy','y-velocity'],['E','TotalEnergy']]:
            self.quan[b]=self.quan[a]
        self.right = p49_eigen.eigen(d, vx, vy, vz, Bx,By,Bz,E, Gamma=Gamma)
        self.speeds = p49_eigen.speeds(d, vx, vy, vz, Bx,By,Bz,E, Gamma=Gamma)
    def perturb(self,base_size=None,pert=1e-6,wave=6,directory=".", write=True):

        self.cubes={}
        map_to_eigen={'d':0,'vx':2,'vy':3,'vz':4,'Bx':-1,'By':5,'Bz':6,'E':1}
        map_to_label ={'d':'density','vx':'x-velocity','vy':'y-velocity','vz':'z-velocity',
                'Bx':'Bx','By':'By','Bz':'Bz','E':'TotalEnergy'}
        face_offset = {'Bx':nar([1,0,0]),'By':nar([0,1,0]),'Bz':nar([0,0,1])}
        for field in ['d','E','vx','vy','vz','Bx','By','Bz']:
            size = base_size+face_offset.get(field,0)
            this_set = np.ones(size)*quan[field]
            if field in ['vx','vy','vz','E']:
                this_set *= self.cubes['density']
            self.cubes[map_to_label[field]]=this_set
        for field in ['d','E','vx','vy','vz','Bx','By','Bz']:
            size = base_size+face_offset.get(field,0)
            amplitude = np.zeros(size)
            amplitude[:size[0]/2,:,:] = 1
            amplitude[size[0]/2:,:,:] = -1
            if field is not 'Bx':
                the_wave = pert*amplitude*self.right[ map_to_eigen[field] ][wave] 
                self.cubes[map_to_label[field]] += the_wave
        for field in ['d','E','vx','vy','vz','Bx','By','Bz']:
            this_filename = "%s/%s_16.h5"%(directory,map_to_label[field])
            if field in ['vx','vy','vz','E']:
                self.cubes[map_to_label[field]] /= self.cubes[map_to_label['d']]
            if write:
                enzo_write.dump_h5(self.cubes[map_to_label[field]],this_filename)

right_fast = p49_eigen.waves(P=0.6)
#right_fast.perturb(base_size=nar([16]*3),pert= 1e-6,wave= 6,directory= '/Users/dcollins/scratch/Paper49b_play/r17d')
print(right_fast['dv'])
"""
below are old things.
"""
def r16_cube_right_fast(directory, write=True):
    #base_size = nar([64]*3)
    base_size = nar([16]*3)
    pert = 1e-6
    wave=6
    quan={}
    d=1.
    vx=vy=vz=0
    Bx, By, Bz = 1.0,1.41421,0.5
    P = 0.6
    Gamma=1.6666666667
    Egas=P/((Gamma-1)*d)
    E = 0.5*(vx*vx+vy*vy+vz*vz)+0.5*(Bx*Bx+By*By+Bz*Bz)/d + Egas
    quan.update({'d':d,'vx':vx,'vy':vy,'vz':vz,'Bx':Bx,'By':By,'Bz':Bz,'P':P,'Gamma':Gamma,'Egas':Egas,'E':E})

    right = p49_eigen.eigen(d, vx, vy, vz, Bx,By,Bz,E, Gamma=Gamma)
    speeds = p49_eigen.speeds(d, vx, vy, vz, Bx,By,Bz,E, Gamma=Gamma)


    cubes={}
    thing={'d':0,'vx':2,'vy':3,'vz':4,'Bx':-1,'By':5,'Bz':6,'E':1}
    thing2={'d':'density','vx':'x-velocity','vy':'y-velocity','vz':'z-velocity',
            'Bx':'Bx','By':'By','Bz':'Bz','E':'TotalEnergy'}
    off = {'Bx':nar([1,0,0]),'By':nar([0,1,0]),'Bz':nar([0,0,1])}
    for field in ['d','E','vx','vy','vz','Bx','By','Bz']:
        size = base_size+off.get(field,0)
        print(field,size)
        this_set = np.ones(size)*quan[field]
        if field in ['vx','vy','vz','E']:
            this_set *= cubes['density']
        cubes[thing2[field]]=this_set
    for field in ['d','E','vx','vy','vz','Bx','By','Bz']:
        size = base_size+off.get(field,0)
        print(field,size)
        amplitude = np.zeros(size)
        amplitude[:size[0]/2,:,:] = 1
        amplitude[size[0]/2:,:,:] = -1
        if field is not 'Bx':
            #this_set += pert*amplitude*right[ thing[field] ][wave]
            the_wave = pert*amplitude*right[ thing[field] ][wave] 
            #this_set = cubes[thing2[field]]+ the_wave
            cubes[thing2[field]] += the_wave
    for field in ['d','E','vx','vy','vz','Bx','By','Bz']:
        this_filename = "%s/%s_16.h5"%(directory,thing2[field])
        if field in ['vx','vy','vz','E']:
            cubes[thing2[field]] /= cubes[thing2['d']]
        if write:
            enzo_write.dump_h5(cubes[thing2[field]],this_filename)

    for a,b in [['d','Density'],['vx','x-velocity'],['vy','y-velocity'],['E','TotalEnergy']]:
        quan[b]=quan[a]
    return cubes, right, quan


#g2, eigen, quan=r16_cube_right_fast('/Users/dcollins/scratch/Paper49b_play/r17c_check', write=True)
#EigenVector=np.zeros([7,7])
#EigenVector[0,6] = 0.8944269506760397 
#EigenVector[1,6] = 4.0249116650296823 
#EigenVector[2,6] = 1.7888514981118755 
#EigenVector[3,6] = -0.8432735801075147 
#EigenVector[4,6] = -0.2981429844604107 
#EigenVector[5,6] = 1.6865448944190704 
#EigenVector[6,6] = 0.5962851678389597 
#for j in range(6):
#    print("E [%d] %0.16f %0.16f"%(j,eigen[j,6], eigen[j,6]-EigenVector[j,6]))
#


def test_eign():
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
