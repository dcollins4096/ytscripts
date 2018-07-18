
from go_lite_pyt3 import *
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
from p49_stuff import *


if 0:
    #Target: in progress
    size = 32
    #k3 = p49_eigen.make_k_freqs(size)
    ks = p49_eigen.make_k_freqs_and_int(size)
    k3 = ks['k_freq']
    kint = ks['k_int']
    k_norm = np.sqrt( k3[0,...]**2 + k3[1,...]**2 + k3[2,...]**2 )
    ok = k_norm > 0
    ampl=np.zeros_like(k_norm)*1j
    theta = np.pi*2*np.random.random(ok.sum())
    phase =  np.exp(theta*1j)
    ampl[ok] = k_norm[ok]**(-5./3)*phase

size = 32
if 1:
    blah=False
    #Target woodshed
    #k3 = p49_eigen.make_k_freqs(size)
    ks = p49_eigen.make_k_freqs_and_int(size)
    k3 = ks['k_freq']
    kint = ks['k_int']
    kint_norm = np.sqrt( kint[0,...]**2 + kint[1,...]**2 + kint[2,...]**2 )
    k_norm = np.sqrt( k3[0,...]**2 + k3[1,...]**2 + k3[2,...]**2 )
    ok = k_norm > 0
    ampl=np.zeros_like(k_norm)*1j
    theta = np.pi*2*np.random.random(ok.sum())
    phase = 1.# np.exp(theta*1j)
    ampl[ok] = k_norm[ok]**(-5./3)*phase

    ampl=np.ones_like(k_norm)*(1+0j)
    this_mask = np.ones_like(ampl,dtype=bool)

    if 0:
        # === target randomly
        this_mask = kint[0,...] == 2
        #this_mask = np.logical_or(this_mask, kint[0,...] ==4)
        #this_mask = np.logical_and(this_mask, kint[1,...] ==2)
        #this_mask = np.logical_and(this_mask, kint[1,...] < 16)
        #this_mask = np.logical_and(this_mask, kint[2,...] == 2)
        #this_mask = k_norm < 5
    if 0:

        # ===  k==2
        # produces f+ f- as expected.  Doesn't produce mixed
        # modes.
        this_mask = np.abs((k_norm*32-2.)) < 0.01/32

    if 0:
        # === k==3, only postive k values.
        # this mixes f+ and f- in a way I don't understand.
        this_mask = np.abs((k_norm*32-3.)) < 0.01/32
        this_mask = np.logical_and(this_mask, kint_norm<4)

    if 0:
        # === Why does this produce f+ and f-?
        #     Still on the kz = 0 plane.
        this_mask[:]=False
        this_mask[0,3,0]=True

        #this_mask[1,1,1]=True
        #this_mask[0,0,3]=True
        #(array([0, 0, 1, 2, 2, 3]),
        # array([0, 3, 2, 1, 2, 0]),
        #  array([3, 0, 2, 2, 1, 0]))

    if 0:
        # === Why does this produce f+ and f-?
        #     Still on the kz = 0 plane.
        this_mask[:]=False
        this_mask[1,2,2]=True

    if 0:
        this_mask = np.logical_and(this_mask,kint[2,...] != 0)

    if 0:
        # === k==3, only postive k values AND non-zero kz.
        # Doesn't mix f+ and f-
        this_mask = np.abs((k_norm*32-3.)) < 0.01/32
        this_mask = np.logical_and(this_mask, kint_norm<4)
        this_mask = np.logical_and(this_mask,kint[2,...] != 0)

    if 0:
        # === all positive k.  Works, no mode mixing.
        this_mask = np.logical_and(this_mask, k3[0,...]>0)
        this_mask = np.logical_and(this_mask, k3[1,...]>0)
        this_mask = np.logical_and(this_mask, k3[2,...]>0)
        
    if 0:
        # === all negative k.  Mixes modes.
        this_mask = np.logical_and(this_mask, k3[0,...]<0)
        this_mask = np.logical_and(this_mask, k3[1,...]<0)
        this_mask = np.logical_and(this_mask, k3[2,...]>0)

    if 1:
        # === all negative k.  Mixes modes.
        this_mask[:] = False
        this_mask[2,2,1]=True


    blank_mask = this_mask == False
    ampl[blank_mask]=0.

    #ampl[1,1,1]=a2[1,1,1]
    #ampl[3,3,3]=a2[3,3,3]
    #ampl[7,8,9]=a2[7,8,9]
if 0:
    import matplotlib.colors as colors
    img_args={ 'origin':'lower','interpolation':'nearest'}
    img_args['norm'] = colors.Normalize(vmin=0,vmax=np.abs(ampl).max())
    p49_plot_tools.six_plot(ampl,oname = 'AMPL.png',img_args=img_args)
    p49_plot_tools.six_plot(k_norm,oname = 'K_NORM.png')
    p49_plot_tools.six_plot(kint_norm,oname = 'K_INT.png')

if 0:
    #works.
    kint = nar([[5,1],[7,0],[8,0]])
    ratio= np.sqrt(2)*(1+1j)
    ampl = nar([1e-6*ratio,0])
    blah=True

if 0:
    #works.
    kint = nar([[2,1],[2,0],[1,0]])
    ratio= 1.#np.sqrt(2)*(1+1j)
    ampl = nar([1e-3*ratio,0])
    blah=True

if  'unlazy' not in dir():
    unlazy = True

if 'ts1' not in dir() or unlazy:
    old_b   = nar([1.0, 1.41421, 0.5])
    WRITE = True
    directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/s01_k53_f+'
    wave = 's+'
    ts1 = p49_eigen.waves(hx=old_b[0],hy=old_b[1],hz=old_b[2],p=0.6,
                          this_wave=wave, form='rb96', HydroMethod=4)
    ts1.rot_write(pert_shape='fft',base_size=nar([size]*3),
                  pert=ampl,directory=directory,
                  k_rot=kint,
                  wave=wave, start=True,write=WRITE, blah=blah)

if unlazy:
    stuff2={'cubes':ts1.cubes,'means':ts1.quan}
    stuff2['ffts']=p49_eigen.get_ffts(stuff2['cubes'],stuff2['means'])
    stuff2['return_system']=p49_eigen.rotate_back(stuff2['ffts'],stuff2['means'])
    stuff2['k_return']=stuff2['return_system'][0]
    stuff2['waves'] = stuff2['return_system'][1].wave_content
    stuff = p49_plot_tools.chomp(directory)

plt.close('all')
if 0:
    p49_plot_tools.print_wave_content(mwut=stuff2['return_system'][1])

analysis={'print_wave':True,
          'plot_fields':True,
          'k_mag':True,
          'k_proj':True}
#analysis={}

od='/Users/dcollins/RESEARCH2/Paper49_EBQU/2018-06-12-p49b/Spectral/'
od += 's01_k53_f+'
if analysis.get('print_wave',False):
    p49_plot_tools.print_wave_content(stuff=stuff)

if analysis.get('plot_fields',False):
    prefix = "%s/%s"%(od,"fields")

    p49_plot_tools.plot_var_fields(stuff=stuff, prefix=prefix, do='slice')

if analysis.get('k_proj',False):
    prefix = "%s/%s"%(od,"K")
    p49_plot_tools.plot_k_proj(stuff=stuff, prefix=prefix)

if analysis.get('k_mag',False):
    prefix = "%s/Mag.png"%od
    p1d=p49_plot_tools.plot_wave_mag(stuff=stuff,output_name=prefix)
    #p49_plot_tools.print_wave_content(stuff=stuff)




if 0:
    for wave in ['f-', 'a-','s-','c','f+','a+','s+']:
        this_fft = stuff2['waves'][wave]
        vvv = np.abs(this_fft).max()
        if vvv > 1e-12:
            max_str = "%5.2e"%vvv
        else:
            max_str = "%5s"%"-"
        print("ampl wave %3s this_max %s "%(wave, max_str))
    

if 0:
    #check the content of the ffts.
    #field_list=['d']
    field_list = ['d','px','py','pz','hx','hy','hz','p']
    for field in field_list:
        if 1:
            print("=== %s xyz ==="%field)
            target_value =  (ampl*ts1.rot[field])[0]
            the_fft = stuff2['ffts'][field]
        if 0:
            print("=== %s FT xyz ==="%field)
            target_value =  (ampl*ts2.rot[field])[0]
            the_fft = ts2.all_hats[field]
        if 0:
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
        kvec = list(zip(k1[0],k1[1],k1[2]))
        #the values
        nozo=nonzero(the_fft)

        #print(sfrm%"+++target+++" + " %s"%cf.format(target_value))
        if 1:
            #number of non-zero k
            print(sfrm%"Nk" + " %d"%nozo.size)
        if 0:
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
                    if False and np.abs(target_value) > 0:
                        print(sfrm%"ratio" + " %s"%cf.format(nozo[n]/target_value))
                    else:
                        print(sfrm%"ratio (zero) "+ "%s"%cf.format(nozo[n]))
