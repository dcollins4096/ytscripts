if 'ef' not in dir():
    execfile('go')
    for i in range(3):
        print("====================")
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
def symmetric(v):
    #for real fft, Ak = A(-k)^*; the negative-phase is the conjugate
    #(that way when you sum, the imaginary parts cancel.)
    #This lets us take an arbitrary K-space signal and ensure it's inverse-fft is real.
    s2=np.zeros_like(v)
    Nx=v.size
    s2[1:Nx/2] = v[1:Nx/2]
    s2[Nx:Nx/2:-1] = v[1:Nx/2].conj()
    s2[0]=v[0]
    return s2

if 1:
    #Rotation test.
    directory = '.'
    name = 'r901'
    wave='f-'
    def stuff(h=None,k_test=None):
        ts = p49_eigen.waves(hx=h[0],hy=h[1],hz=h[2],p=0.6,this_wave=wave) #, form='rb96')
        #k_test = nar([[1.,0.],[0.,0.],[0.,1]])
        ampl = nar([1e-6,0])
        kint = k_test.astype(int)
        ts.rotate(kint)
        def ugprint(label,thing, form='%5.2f'):
            stuff1 = "  %s"%form
            stuffN = stuff1*len(thing)

            print(label+stuffN%tuple(thing))
        print(" ")
        ugprint('h_unit ',ts.h_unit)
        print(" ")
        ugprint('a_unit0',ts.a_unit[:,0])
        ugprint('b_unit0',ts.b_unit[:,0])
        ugprint('c_unit0',ts.c_unit[:,0])
        ugprint('B0hat0 ',ts.B0hat[:,0])
        ugprint('dh     ',(ts.rot['hx'][0],ts.rot['hy'][0],ts.rot['hz'][0]))
        ugprint('dv     ',(ts.rot['vx'][0],ts.rot['vy'][0],ts.rot['vz'][0]))
        ugprint('drho,de',(ts.rot['d'][0],ts.rot['e'][0]))

        #ugprint('a_unit1',ts.a_unit[:,1])
        #ugprint('b_unit1',ts.b_unit[:,1])
        #ugprint('c_unit1',ts.c_unit[:,1])
        #ugprint('B0hat1 ',ts.B0hat[:,1])
        #ugprint('db     ',(ts.rot['hx'][1],ts.rot['by'][1],ts.rot['bz'][1]))
        #ugprint('vb     ',(ts.rot['vx'][1],ts.rot['vy'][1],ts.rot['vz'][1]))
        return ts
    if 0:
        print("Test 1: H in xy plane.  full \pi/2 rotations to the coordinate axes.")
        print("        Magnetic vectors should cyclicly permute.")
        t1 =stuff(h=[1.0,0.5,0.0],k_test=nar([[1.,0.],[0.,0.],[0.,1]]))
        t1a=stuff(h=[0.0,1.0,0.5],k_test=nar([[0.,0.],[1.,0.],[0.,1]]))
        t1b=stuff([0.5,0.0,1.0],nar([[0.,0.],[0.,0.],[1.,1]]))
    if 0:
        print("Test 2: H in xz plane. full \pi/2 rotations to the coordinate axes.")
        print("        Magnetic vectors should cyclicly permute.")
        t2 =stuff(h=[1.0,0.0,0.5],k_test=nar([[1.,0.],[0.,0.],[0.,1]]))
        t2a=stuff(h=[0.5,1.0,0.0],k_test=nar([[0.,0.],[1.,0.],[0.,1]]))
        t2b=stuff(h=[0.0,0.5,1.0],k_test=nar([[0.,0.],[0.,0.],[1.,1]]))
    if 0:
        print("Test 3: H in xyz plane. full \pi/2 rotations to the coordinate axes.")
        print("        b_unit")
        t3 =stuff(h=[1.0,0.5,0.2],k_test=nar([[1.,0.],[0.,0.],[0.,1]]))
        t3a=stuff(h=[0.2,1.0,0.5],k_test=nar([[0.,0.],[1.,0.],[0.,1]]))
        t3b=stuff(h=[0.5,0.2,1.0],k_test=nar([[0.,0.],[0.,0.],[1.,1]]))
    if 1:
        print("Test 4: our standard test")
        t3 =stuff(h=[1.0,0.5,0.2],k_test=nar([[1.,0.],[0.,0.],[0.,1]]))
        t3a=stuff(h=[0.2,1.0,0.5],k_test=nar([[0.,0.],[1.,0.],[0.,1]]))
        t3b=stuff(h=[0.5,0.2,1.0],k_test=nar([[0.,0.],[0.,0.],[1.,1]]))

   #print("h_unit",ts.h_unit)
   #print("a_unit",ts.a_unit)
   #print("b_unit",ts.b_unit)
   #print("c_unit",ts.c_unit)
   #print("B0hat",ts.B0hat)

if 0:
    this_system = p49_eigen.waves(bx=1.0,by=0.0,bz=0.0,p=0.6,this_wave='f-')
    k_test = nar([[1.,1.],[1.,0.],[0.,1]])
    #ampl = nar([1e-2,1e-2])
    ampl = nar([1e-2,0])
    #this_system.rotate(k_test)
    #kint = k_test.astype(int)
    this_system.rot_write(k_test,base_size=nar([16]*3),pert=ampl,directory='/Users/dcollins/scratch/Paper49b_play/r103a_f-_110')
if 0:
    this_system = p49_eigen.waves(bx=1.0,by=0.0,bz=0.0,p=0.6,this_wave='a+')
    k_test = nar([[2.,1.],[1.,1.],[0.,1]])
    ampl = nar([1e-6,1e-6])
    kint = k_test.astype(int)
    this_system.rot_write(k_test,base_size=nar([16]*3),pert=ampl,directory='/Users/dcollins/scratch/Paper49b_play/r101_two')
if 0:
    this_system = p49_eigen.waves(bx=1.0,by=0.0,bz=0.0,p=0.6,this_wave='a+')
    k_test = nar([[1.,1.],[1.,0.],[0.,1]])
    #ampl = nar([1e-2,1e-2])
    ampl = nar([0,1e-2])
    kint = k_test.astype(int)
    this_system.rot_write(k_test,base_size=nar([16]*3),pert=ampl,directory='/Users/dcollins/scratch/Paper49b_play/r102b_wave_101')

if 0:
    this_system = p49_eigen.waves(bx=1.0,by=0.0,bz=0.0,p=0.6,this_wave='a+')
    k_test = nar([[1.,1.],[1.,0.],[0.,1]])
    ampl = nar([1e-2,1e-2])
    ampl = nar([0,1e-2])
    kint = k_test.astype(int)
    this_system.rot_write(k_test,base_size=nar([16]*3),pert=ampl, write=False) #,directory='/Users/dcollins/scratch/Paper49b_play/r104_linear_combo')
    ampl = nar([1e-2,0])
    this_system.rot_write(k_test,base_size=nar([16]*3),pert=ampl,directory='/Users/dcollins/scratch/Paper49b_play/r104_linear_combo', wave='a-')

if 0:
    this_system = p49_eigen.waves(bx=1.0,by=1.41421,bz=0.5,p=0.6,this_wave='f-') #, form='rb96')
    k_test = nar([[1.,1.],[0.,0.],[0.,1]])
    ampl = nar([1e-2,0])
    ampl = nar([1e-6,0])
    #ampl = nar([0,1e-2])
    kint = k_test.astype(int)
    directory = '/Users/dcollins/scratch/Paper49b_play/r201_rb96_f-'
    name = 'r201'
    this_system.rot_write(k_test,base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave='f-')
    temp_taxi = taxi.taxi(directory=directory, name = name, frames=[0,50],axis=[0,1,2])
    field_list = ['Density','x-velocity','y-velocity','z-velocity','Bx','By','Bz','pressure','TotalEnergy']
    temp_taxi.Colorbar='fixed'
    temp_taxi.operation='CenterSlice'
    temp_taxi.fields=field_list
    for field in field_list:
        ifield = field
        if ifield in ['Bx','By','Bz']:
            ifield = field.lower()
        if ifield == 'pressure':
            ifield='p'
        temp_taxi.slice_zlim[field] =  [this_system.quan[ifield]-2e-6,this_system.quan[ifield]+2e-6]
        temp_taxi.set_log[field]=False


