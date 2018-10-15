
#if 'ef' not in dir():
#    execfile('go')
#    for i in range(3):
#        print("====================")
from go_lite_pyt3 import *
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)

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
    return None
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
    this_system = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, HydroMethod = 6) #, form='rb96')
    ampl = nar([1e-6,0])
    this_system.rot_write(pert_shape='square_x',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave)
    #r401=eigen_taxi(this_system,name)
    
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
    ts = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96',HydroMethod=6)
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
if 1:
    #r701b Hydro 4
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r701b_hydro4'
    name = 'r701'
    wave='f-'
    ts = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96',HydroMethod=4)
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
    ts801 = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
    k_test = nar([[1.,1.],[1.,0.],[0.,1]])
    ratio = 1.0# ts.speeds['cf']/ts.speeds['aa']
    ampl = nar([1e-6,0])
    ts801.rot_write(pert_shape='fft',base_size=nar([16]*3),pert=ampl,directory=directory,
                          wave=wave,k_rot=k_test)
    r801=eigen_taxi(ts801,name)

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

if 0:
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
if 0:
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

if 0:
    #IN PROCESS
    #rA01
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/rA01_rb96_110_f-'

if 0:
    #rB01 rB01_rb_several
    #messing around with rotated vectors.
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/rB01_rb_several'
    name = 'rB01'
    wave='s-'
    #k_test = nar([[1.,0.],[0.,0.],[0.,1]])
    k_test = nar([[3.,0.],[2.,0.],[0.,1]])
    ratio = 1.0# ts.speeds['cf']/ts.speeds['aa']
    size=32
    ampl = nar([1e-6*ratio,0])
    def rot(a,b,theta):
        ahat = np.cos(theta)*a-np.sin(theta)*b
        bhat = np.sin(theta)*a+np.cos(theta)*b
        return ahat,bhat
    old_b   = rot(1.0, 1.41421, 0)
    old_bz = 0.5
    old_K   = rot(5, 0, 0)
    old_K   = rot(5, 8, 0)
    k1 = nar([[old_K[0],0.],[old_K[1],0.],[0.,1.]])
#   ts1 = p49_eigen.waves(hx=old_b[0],hy=old_b[1],hz=old_bz,p=0.6,this_wave=wave, form='rb96', HydroMethod=4)
#   ts1.rot_write(pert_shape='fft',base_size=nar([32]*3),pert=ampl,directory=directory,
#                        wave='f-',k_rot=k1,
#                            start=True,write=False)
    this_theta = 0.0 #np.arctan(3./4)
    new_b   = rot(1.0, 1.41421, this_theta)
    new_K   = rot(-1, -1, this_theta)
    new_bz = old_bz
    k2 = nar([[new_K[0],0.],[new_K[1],0.],[0.,1.]])
    k2 = nar([[5,1],[7,0],[8,0]])
    ts2 = p49_eigen.waves(hx=new_b[0],hy=new_b[1],hz=new_bz,p=0.6,this_wave=wave, form='rb96', HydroMethod=4)
    directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/s01_k53_f+'
    ts2.rot_write(pert_shape='fft',base_size=nar([size]*3),pert=ampl,directory=directory,
                         wave='s-',k_rot=k2,
                             start=True,write=True)
    #ts1.rot_write(pert_shape='fft',base_size=nar([32]*3),pert=ampl,directory=directory,
    #                     wave='s-',k_rot=k2,
    #                         start=False,write=False)
    #ts1.rot_write(pert_shape='fft',base_size=nar([32]*3),pert=ampl,directory=directory,
    #                      wave='f-',
    #                      k_rot=nar([[3,0],[2,0],[0,1]]),start=True,write=False)
#    stuff1={'cubes':ts1.temp_cubes,'means':ts1.quan}
#    stuff1['ffts']=p49_eigen.get_ffts(stuff1['cubes'])
    stuff2={'cubes':ts2.temp_cubes,'means':ts2.quan}
    stuff2['ffts']=p49_eigen.get_ffts(stuff2['cubes'],stuff2['means'])

if 0:
    #they rotate back properly.
    ts1.project_to_waves(k1[:,0],ts1.pull_rot(0),means={})
    ts2.project_to_waves(k2[:,0],ts2.pull_rot(0),means={})
    def slim(v):
        if v > 1e-12:
            return "%10.2e"%v
        else:
            return "%10s"%"-"
    for wave in ['f-', 'a-','s-','c','f+','a+','s+']:
        print("=== W %3s %s %s"%(wave,slim(ts1.wave_content[wave]), slim(ts2.wave_content[wave])))


if 0:
    field_list = ['d','vx','vy','vz','hx','hy','hz','p']
    kall,mwut=p49_eigen.rotate_back(stuff2['ffts'],stuff2['means'])
    sfrm = "%12s"
    this_ts = ts2
    these_means = stuff2['means']

if 0:
    for field in field_list:
        for wave in ['s+','s-']:
            this_wave = mwut.wave_content[wave]
            k1 = nz(this_wave)
            nozo=nonzero(this_wave)
            kvec = zip(k1[0],k1[1],k1[2])
            print("wave %3s f %3s %s"%(wave,field,ampl_str(mwut.hat_system.left[wave][field][k1])))
            print(" wave_frame    %s"%(ampl_str(0.5e6*mwut.wave_frame[field][k1])))
            print(" k             %s\n"%(zip(kall[0,...][k1],kall[1,...][k1],kall[2,...][k1])))

if 0:
    for wave in ['f-', 'a-','s-','c','f+','a+','s+']:
        this_wave = mwut.wave_content[wave]
        k1 = nz(this_wave)
        nozo=nonzero(this_wave)
        kvec = zip(k1[0],k1[1],k1[2])
        print("=== %s ==="%wave)
        for n in range(len(kvec)):
            print("  "+str(kvec[n]))
            print("  "+cf.format(nozo[n]))

        #print("%5s %s"%(wave,kvec))
        #print("%5s %s"%("",str(nozo)))

if 0:
    #check the content of the ffts.
    #field_list=['d']
    for field in field_list:
        if 0:
            print("=== %s xyz ==="%field)
            target_value =  (ampl*ts2.rot[field])[0]
            the_fft = stuff2['ffts'][field]
        if 0:
            print("=== %s FT xyz ==="%field)
            target_value =  (ampl*ts2.rot[field])[0]
            the_fft = ts2.all_hats[field]
        if 1:
            print("=== %s FT^-1(FT) ==="%field)
            target_value =  (ampl*ts2.rot[field])[0]
            the_fft = ts2.temp_right_back[field]
        if 0:
            print("=== %s abc ==="%field)
            #target_value = ts2.hat_system.quan[field]
            target_value = ts2.hat_system.right['s-'][field][0]
            the_fft = mwut.wave_frame[field]
        #the vectors
        k1 = nz(the_fft)
        kvec = zip(k1[0],k1[1],k1[2])
        #the values
        nozo=nonzero(the_fft)

        print(sfrm%"+++target+++" + " %s"%cf.format(target_value))
        if 0:
            #number of non-zero k
            print(sfrm%"Nk" + " %d"%nozo.size)
        if 1:
            print(sfrm%"size " + str(the_fft.shape))
            for n ,k in enumerate(kvec):
                print(sfrm%"k  "+str(k))
        if 0:
            #diagnostics about non-zero k, values, and expected values
            for n ,k in enumerate(kvec):
                is_mean = False
                if (nar(k) == 0 ).all():
                    is_mean = True
                #What vectors are populated?  
                #        [0,0,0] if there's a mean.
                #                others for the actual wave
                print(sfrm%"k vector"+ " %s"%str(k))
                    
                #what is the target?
                #print(sfrm%"eigen" + " %s"%ampl_str([target_value]))


                if is_mean:
                    print(sfrm%"quan" +" %s"%cf.format(this_ts.quan[field]))
                    print(sfrm%"mean" +" %s"%cf.format(these_means[field]))
                    print(sfrm%"amplitude" +" %s"%cf.format(nozo[n]))
                else:
                    #What's the amplitude of the values?
                    #        Mean, if there is one.
                    #        Matches the target
                    print(sfrm%"amplitude" +" %s"%cf.format(nozo[n]))

                    #How did we do?
                    if np.abs(target_value) > 0:
                        print(sfrm%"ratio" + " %s"%cf.format(nozo[n]/target_value))
                    else:
                        print(sfrm%"ratio (zero) "+ "%s"%cf.format(nozo[n]))
if 0:
    def maxer(arr):
        out = "-"
        themax = np.abs(arr).max() 
        if themax> 1e-12:
            out = "%5.2e"%themax
        return out

    print("===== wave content =====")
    for wave in ['f-', 'a-','s-','c','f+','a+','s+']:
        print( "%5s %5s"%(wave,maxer(mwut.wave_content[wave])))


if 0:
    these_means = stuff2['means']
    these_ffts  = p49_eigen.get_ffts(stuff2['cubes'], these_means)
    #kall,lwut=p49_eigen.rotate_back(these_ffts, these_means)
    kall,lwut=p49_eigen.rotate_back(these_ffts, stuff2['means']) #_means)
    kall,mwut=p49_eigen.rotate_back(stuff2['ffts'],stuff2['means'])
    for wave in ['f-', 'a-','s-','c','f+','a+','s+']:
        print( "%5s %5s"%(wave,maxer(lwut.wave_content[wave])))

    print('aaa')
    p49_plot_tools.print_wave_content(mwut=lwut)
    print('bbb')
    p49_plot_tools.print_wave_content(mwut=mwut)

if 0:
    #The content is not right; too much in t2 in the 'wrong' mode.
    reload(p49_plot_tools)
    #p49_plot_tools.print_wave_content(stuff=stuff1)
    p49_plot_tools.print_wave_content(stuff=stuff2)


if 0:
    import p49_plot_tools
    reload(p49_plot_tools)
    plot_directory = "/Users/dcollins/RESEARCH2/Paper49_EBQU/2018-06-12-p49b/multi_wave/"
    field_list = ['d','vx','vy','vz','hx','hy','hz','p']
    for field in field_list:
        print("== %s =="%field)
        print("   1r rot %5.2e max ft %0.2e"%(ts1.rot[field][0], np.abs(ts1.all_p[field]).max()))
        print("   2r rot %5.2e max ft %0.2e"%(ts2.rot[field][0], np.abs(ts2.all_p[field]).max()))
if 0:
    p49_plot_tools.print_wave_content(stuff=stuff1)
    p49_plot_tools.print_wave_content(stuff=stuff2)
    prefix="rB01_0661"
    p49_plot_tools.plot_wave_mag(stuff=stuff1,output_name=plot_directory+prefix+"_Kmag.png")
    p49_plot_tools.plot_k_proj(stuff=stuff1,prefix="%s%s"%(plot_directory,prefix))
    prefix="rB01_0662"
    p49_plot_tools.plot_wave_mag(stuff=stuff2,output_name=plot_directory+prefix+"_Kmag.png")
    p49_plot_tools.plot_k_proj(stuff=stuff2,prefix="%s%s"%(plot_directory,prefix))

    #ts2 = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96', HydroMethod=4)
    #ts2.rot_write(pert_shape='fft',base_size=nar([32]*3),pert=ampl,directory=directory,
    #                     wave=wave,k_rot=k_test,start=True)
    #ts2.rot_write(pert_shape='fft',base_size=nar([32]*3),pert=ampl,directory=directory,
    #                      wave='s+',k_rot=nar([[3,4],[2,4],[0,1]]),start=False)
    #ts3 = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96', HydroMethod=4)
    #ts3.rot_write(pert_shape='fft',base_size=nar([32]*3),pert=ampl,directory=directory,
    #                      wave=wave,k_rot=k_test,start=True)
    #ts3.rot_write(pert_shape='fft',base_size=nar([32]*3),pert=ampl,directory=directory,
    #                      wave='s+',k_rot=nar([[3,4],[2,4],[0,1]]),start=True)
    #print('bigger.')
    #rB01=eigen_taxi(ts,name)
