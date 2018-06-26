
if 'ef' not in dir():
    execfile('go')
    for i in range(3):
        print("====================")
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)


def prnt1(thing):
    for field in ['vx','vy','vz','hx','hy','hz']:
        print("%4s %5.2f"%(field,thing[field]))
def prnt2(thing):
    for field in ['vx','vy','vz','hx','hy','hz']:
        print("%4s %5.2f %5.2f"%(field,thing[field][0],thing[field][1]))
def prnt1b(thing):
    for field in ['d','p','vx','vy','vz','hx','hy','hz']:
        print("%4s %5.2f"%(field,thing[field]))
def prnt2b(thing):
    for field in ['d','p','vx','vy','vz','hx','hy','hz']:
        print("%4s %5.2f %5.2f"%(field,thing[field][0],thing[field][1]))

reload(p49_eigen)

if 0:
    #r701 r701_rb_fft_f-
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r701_rb96_fft_f-'
    name = 'r701'
    wave='f-'
    ts701 = p49_eigen.waves(hx=1.0,hy=1.0,hz=1.0,p=0.6,this_wave=wave, form='rb96')
    k_test = nar([[1.,0.],[0.,0.],[0.,1]])
    ratio = ts701.speeds['cf']/ts701.speeds['aa']
    ampl = nar([1e-6*ratio,2.3])
    ampl = nar([1,0.2])
    ampl = nar([1,0.])
    ts701.rot_write(pert_shape='fft',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test,write=False)
    ret701 = p49_eigen.waves(form=ts701.form,**ts701.quan)

if 0:
    #eigen vector check
    ret701.check_orthonormality(k_test,ts701.rot)
if 0:
    #rotation check; dummy is dumb, rotten rotates, stillrotten is dummy back again.
    dummy={'vx':0.1,'vy':0.2,'vz':0.3,'hx':0,'hy':0,'hz':0}
    rotten = ts701.rotate_to_xyz(dummy)
    stillrotten = ts701.rotate_to_abc(rotten)
    print("=== original ===")
    prnt1(dummy)
    print("=== two rotated ===")
    prnt2(rotten)
    print("=== rotated back ===")
    prnt2(stillrotten)

if 0:
    #rot is the eigen vectors along k_test, rotated to physics.
    #Successfully extracts just the f- wave.
    ret701 = p49_eigen.waves(form=ts701.form,**ts701.quan)
    ret701.fields_to_wave_frame(k_test,ts701.rot)
    
    print("=== eigen: truth ===")
    prnt1b(ts701.right[wave])
    print("=== real space: along k ===")
    prnt2b(ts701.rot)
    print("=== wave space: these should match ===")
    prnt2b(ret701.wave_frame)


    #ret701.project_to_waves(k_test, ts701.rot)
    zero = {}
    for frame in ts701.rot:
        zero[frame] = 0.
    ret701.project_to_waves(k_test, ts701.rot, means=zero)
    print("=== wave content: ===")
    for wave in ret701.wave_content:
        print("%5s "%wave, ret701.wave_content[wave])
    # 1.) start with rotated field; t701.rot, real space.
    # 2.) given k, rotate back, project.

if 1:
    wave='f-'
    tsfft = p49_eigen.waves(hx=1.0,hy=1.0,hz=1.0,p=0.6,this_wave=wave, form='rb96')
    k_test = nar([[1.,0.],[0.,0.],[0.,1]])
    ampl = nar([1,1.])
    nn=5
    tsfft.rot_write(pert_shape='fft',base_size=nar([nn]*3),pert=ampl,directory='',
                          wave=wave,k_rot=k_test,write=False)
    k_test = nar([[1.,1.],[0.,0.],[0.,1]])
    tsfft.rot_write(pert_shape='fft',base_size=nar([nn]*3),pert=ampl,directory='',
                          wave='a+',k_rot=k_test,write=False)

if 1:
    field_list = ['d','vx','vy','vz','hx','hy','hz','p']
    def nz(field):
        nz = np.abs(field) > 1e-13
        return field[nz]


    if 0:
        #target: not working
        these_means = stuff['cubes']
        these_ffts = ffts
    if 0:
        #really stupid test.
        field_list = ['d','vx','vy','vz','hx','hy','hz','p']
        these_cubes={}
        for field in field_list:
            these_cubes
        fobase_k = np.zeros([4,4,4])*1j
    if 0:
        #great test
        #test one: actual FFTs and cubes
        these_cubes={}
        more_fft={}
        for a,b in [['d','density'],['vx','x-velocity'],['hx','Bx'],['hy','By'],['hz','Bz'],
                ['vz','z-velocity'],['vy','y-velocity'], ['p','GasPressure']]:
            these_cubes[a] = tsfft.cubes[b][:nn,:nn,:nn]#- tsfft.quan[a]

        #get_ffts subtracts mean, converts to conserved space.
        these_ffts = p49_eigen.get_ffts(these_cubes, tsfft.quan)
    if 0:
        #new thing, in preocess.
        directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r701_rb96_fft_f-'
        frame = 0
        ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
        stuff = p49_eigen.get_cubes_cg(ds)

    if 0:
        ffts = p49_eigen.get_ffts(stuff['cubes'])

    print_fields = False
    print_waves = True
    kall,wut=p49_eigen.rotate_back(these_ffts, these_means)
    fl =  np.zeros_like(wut.wave_frame['d']).astype('bool')
    if print_fields:
        for field in wut.wave_frame:
            print(" ===== %s ===="%field)
            thisthing =  wut.wave_frame[field]
            thisthing =  wut.dumb[field]
            this_bool = np.abs(thisthing) > 1e-13
            #fl = np.logical_or(fl, this_bool)
            nonzeros = len( this_bool )
            print("  eigen    %s"%str(tsfft.right['f-'][field]))
            print("  rot      %s"%str(tsfft.rot[field]))
            print("all_hat  %3s %s"%(field, nz(tsfft.all_hats[field])))
            aaa = these_ffts[field] #is good.
            print("also fft input  k %3s %s"%(field, str(nz(aaa).size)))
            print("this wave frame k %3s %s"%(field, str(nz(thisthing).size)))
    if print_waves:
        for wave in wut.wave_content:
            thisthing = wut.wave_content[wave]
            bang_or_not = ""
            if ( np.abs(thisthing)>1e-13).sum() > 0:
                bang_or_not = "!!!"*8
            print("=== Wave %s %s"%(wave, bang_or_not))
            s1 = str(nz(thisthing.real).size)
            s2 = str(nz(thisthing.imag).size)
            print("wave real nz %s imag %s"%(s1,s2))
        

#
# live test.
# 

