
if 'ef' not in dir():
    execfile('go')
    for i in range(3):
        print("====================")
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)

def uniform(ds=None,cg=None,axis=0,fields=['By'] ):
    if cg is None:
        cg = ds.covering_grid(0,[0.0]*3,[16]*3,fields=fields)
    sl_zero = [slice(None)]*3
    sl_zero[axis] = slice(0,1)
    total_total = 0.0
    for field in fields:
        this_cube = cg[field].v
        this_std = np.std( this_cube,axis=axis)
        full_std = np.std( this_cube)
        total_std = np.sum(np.abs(this_std))
        total_total += total_std
        print("Total std, %13s, %0.3e, full %0.3e"%(field, total_std, full_std))
        #mean = np.mean( np.abs(this_cube), axis=axis)
        #print( this_cube[field][sl_zero])
    if len(fields) > 1:
        print("Total std, all, %0.3e"%(total_total))

    return cg, this_std, total_std

    
print("\ncg, this_std,total = uniform(r701.ds, axis=1, fields=['density','Bx','By','Bz','TotalEnergy','x-velocity','y-velocity','z-velocity'])\n")
def eigen_taxi(this_system, name='p49b_eigen_car', directory=None):
    if directory is None:
        directory = this_system.directory
    temp_taxi = taxi.taxi(directory=directory, name = name, frames=[0,50],axis=[0,1,2])
    field_list = ['Density','x-velocity','y-velocity','z-velocity','Bx','By','Bz','pressure','TotalEnergy']
    temp_taxi.fields=field_list
    temp_taxi.Colorbar='fixed'
    temp_taxi.operation='CenterSlice'
    for field in field_list:
        ifield = {'Bx':'hx','By':'hy','Bz':'hz','pressure':'p'}.get(field,field)
        temp_taxi.slice_zlim[field] =  [this_system.quan[ifield]-2e-6,this_system.quan[ifield]+2e-6]
        temp_taxi.set_log[field]=False
    return temp_taxi

if 0:
    #r401 r401_rj95_sq_f-
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r401_rj95_sq_f-'
    name = 'r401'
    wave='f-'
    this_system = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave) #, form='rb96')
    ampl = nar([1e-6,0])
    this_system.rot_write(pert_shape='square_x',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave)
    r401=eigen_taxi(this_system,name)
    
if 0:
    #r501 r501_rj95_fft_f-
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r501_rj95_fft_f-'
    name = 'r501'
    wave='f-'
    ts501 = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave) #, form='rb96')
    k_test = nar([[1.,1.],[0.,0.],[0.,1]])
    ampl = nar([1e-6,0])
    ts501.rot_write(pert_shape='fft',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test )
    r501=eigen_taxi(ts501,name)

if 0:
    #r601 r601_rb_sq_f-
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r601_rb96_sq_f-'
    name = 'r601'
    wave='f-'
    ts = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
    k_test = nar([[1.,1.],[0.,0.],[0.,1]])
    ratio = ts.speeds['cf']/ts.speeds['aa']
    ampl = nar([1e-6*ratio,0])
    ts.rot_write(pert_shape='square_x',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test)
    r601=eigen_taxi(ts,name)
if 0:
    #r701 r701_rb_fft_f-
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r701_rb96_fft_f-'
    name = 'r701'
    wave='f-'
    ts = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
    k_test = nar([[1.,1.],[0.,0.],[0.,1]])
    ratio = ts.speeds['cf']/ts.speeds['aa']
    ampl = nar([1e-6*ratio,0])
    ts.rot_write(pert_shape='fft',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test)
    r701=eigen_taxi(ts,name)

if 0:
    #NOT WORKING
    #r801 r801_rj_110_f-
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r801_rj95_110_f-'
    name = 'r801'
    wave='f-'
    ts = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave)#, form='rb96')
    k_test = nar([[1.,1.],[1.,0.],[0.,1]])
    ratio = 1.0# ts.speeds['cf']/ts.speeds['aa']
    ampl = nar([1e-6,0])
    ts.rot_write(pert_shape='fft',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test)
    r801=eigen_taxi(ts,name)

if 0:
    #r901 r901_rj_010_f-
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r901_rj95_010_f-'
    name = 'r901'
    wave='f-'
    ts901 = p49_eigen.waves(hx=0.5,hy=1.0,hz=1.41421,p=0.6,this_wave=wave)#, form='rb96')
    k_test = nar([[0.,1.],[1.,0.],[0.,1]])
    ratio = 1.0# ts.speeds['cf']/ts.speeds['aa']
    ampl = nar([1e-6,0])
    #ts901.rotate(k_test)
    ts901.rot_write(pert_shape='fft',base_size=nar([16]*3),pert=ampl,directory=directory, wave=wave,k_rot=k_test)
    #r901=eigen_taxi(ts901,name)

if 0:
    #r402 r402_rj95_sq_f+
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r402_rj95_sq_f+'
    name = 'r402'
    wave='f+'
    this_system = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave) #, form='rb96')
    ampl = nar([1e-6,0])
    this_system.rot_write(pert_shape='square_x',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave)
    r402=eigen_taxi(this_system,name)

if 1:
    #r602 r602_rb_sq_f+
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r602_rb96_sq_f+'
    name = 'r602'
    wave='f+'
    ts = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
    k_test = nar([[1.,1.],[0.,0.],[0.,1]])
    ratio = ts.speeds['cf']/ts.speeds['aa']
    ampl = nar([1e-6*ratio,0])
    ts.rot_write(pert_shape='square_x',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test)
    r602=eigen_taxi(ts,name)

if 0:
    #r603 r603_rb_sq_s-
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r603_rb96_sq_s-'
    name = 'r603'
    wave='s-'
    ts603 = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
    k_test = nar([[1.,1.],[0.,0.],[0.,1]])
    ratio = 1.0#  ts603.speeds['cf']/ts603.speeds['aa']
    ampl = nar([1e-6*ratio,0])
    ts603.rot_write(pert_shape='square_x',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test)
    r603=eigen_taxi(ts603,name)

if 0:
    #r604 r603_rb_sq_s+
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r604_rb96_sq_s+'
    name = 'r604'
    wave='s+'
    ts = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
    k_test = nar([[1.,1.],[0.,0.],[0.,1]])
    ratio = 1.0# ts.speeds['cf']/ts.speeds['aa']
    ampl = nar([1e-6*ratio,0])
    ts.rot_write(pert_shape='square_x',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test)
    r604=eigen_taxi(ts,name)
if 0:
    #r605 r605_rb_sq_a-
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r605_rb96_sq_a-'
    name = 'r605'
    wave='a-'
    ts = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
    k_test = nar([[1.,1.],[0.,0.],[0.,1]])
    ratio = 1.0# ts.speeds['cf']/ts.speeds['aa']
    ampl = nar([1e-6*ratio,0])
    ts.rot_write(pert_shape='square_x',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test)
    r605=eigen_taxi(ts,name)
if 1:
    #r606 r606_rb_sq_a+
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r606_rb96_sq_a+'
    name = 'r606'
    wave='a+'
    ts = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
    k_test = nar([[1.,1.],[0.,0.],[0.,1]])
    ratio = 1.0# ts.speeds['cf']/ts.speeds['aa']
    ampl = nar([1e-6*ratio,0])
    ts.rot_write(pert_shape='square_x',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test)
    r605=eigen_taxi(ts,name)
