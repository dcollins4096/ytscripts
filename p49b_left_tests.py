
if 'ef' not in dir():
    execfile('go')
    for i in range(3):
        print("====================")
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)



reload(p49_eigen)

if 0:
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r701_rb96_fft_f-'
    frame = 0
    ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
    stuff = p49_eigen.get_cubes_cg(ds)

if 0:
    ffts = p49_eigen.get_ffts(stuff['cubes'])

if 0:
    whatever=p49_eigen.rotate_back(ffts,stuff['means'])

if 1:
    #r701 r701_rb_fft_f-
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r701_rb96_fft_f-'
    name = 'r701'
    wave='f-'
    ts701 = p49_eigen.waves(hx=1.0,hy=1.0,hz=1.0,p=0.6,this_wave=wave, form='rb96')
    k_test = nar([[1.,0.],[0.,0.],[0.,1]])
    ratio = ts701.speeds['cf']/ts701.speeds['aa']
    ampl = nar([1e-6*ratio,2.3])
    ampl = nar([1,0.2])
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
if 1:
    ret701 = p49_eigen.waves(form=ts701.form,**ts701.quan)
    ret701.fields_to_wave(k_test,ts701.rot)

if 0:
    #more complete test
    perturbed = {}
    zeros = {}
    for field in ts701.rot:
        #perturbed[field] = np.zeros(2)+ts701.quan[field] + np.ones(2)*0.1*ts701.right['f+'][field] #ampl*ts701.rot[field]
        perturbed[field] = ampl*ts701.rot[field]

    field_list = ['d','vx','vy','vz','hx','hy','hz','p']
    means_zero={}
    for field in field_list:
        means_zero[field]=np.zeros(2)
    ret701.project_to_waves(k_test,perturbed, means=ts701.quan)
    for wave in ret701.wave_content: #['f-', 'a-','s-','c','f+','a+','s+']:
        print( "%5s "%wave, ret701.wave_content[wave])
    #print(ret701.wave_content)


