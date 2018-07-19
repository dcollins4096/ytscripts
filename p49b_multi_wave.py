
from go_lite_pyt3 import *
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
from p49_print_tools import *
plt.close('all')

#size = 128
size = 32
sim = 'rx11'
#sim='rx06_like'

blah=False
if 1:
    #setup FFT init
    blah=False
    #Target woodshed
    #k3 = p49_eigen.make_k_freqs(size)
    ks = p49_eigen.make_k_freqs_and_int(size)
    k3 = ks['k_freq']
    kint = ks['k_int']
    kint_norm = np.sqrt( kint[0,...]**2 + kint[1,...]**2 + kint[2,...]**2 )
    k_norm = np.sqrt( k3[0,...]**2 + k3[1,...]**2 + k3[2,...]**2 )
    ok = np.logical_and(kint[2,...]>0,k_norm > 0)


old_b   = nar([1.0, 1.41421, 0.5])
if sim=='rx11':
    #rx10: ensemble of waves, random phase.
    wavelist = ['f-','s+']
    amplitudes={}
    for wave in wavelist:
        scale=1e-3
        amplitudes[wave]=np.zeros_like(k_norm)*(1+0j)
        #theta = np.pi*2*np.random.random(ampl.shape)
        #phase = np.exp(theta*1j)

        #ampl[mask] = k_norm[mask]**(-5./3)*phase[mask]
        #ampl[mask] *= 1e-3/ampl[mask].max()
    if 'f-' in wavelist:
        mask=ok
        mask[:]=False
        mask[1,2,3]=True
        amplitudes['f-'][mask]=scale
    if 's+' in wavelist:
        mask=ok
        mask[:]=False
        mask[2,1,2]=True
        amplitudes['s+'][mask]=scale
if sim=='rx06_like':
    amplitude=1e-3
    ampl=np.zeros_like(k_norm)*(1+0j)
    old_b   = nar([1.0, 0.5, 1.41421])
    ampl[2,1,2]=amplitude
    amplitudes={'a-':ampl}
out_dir_b= '/Users/dcollins/scratch/Paper49b_play/Spectral/rx11_mixed_1'
out_dir_s= '/Users/dcollins/scratch/Paper49b_play/Spectral/rx11_mixed_s+'
out_dir_f= '/Users/dcollins/scratch/Paper49b_play/Spectral/rx11_mixed_f-'

if 'unlazy' not in dir():
    unlazy=True
if 'tsb' not in dir() or unlazy:
    tsb = p49_eigen.waves(hx=old_b[0],hy=old_b[1],hz=old_b[2],p=0.6,
                           form='rb96', HydroMethod=4)
    blah=False
    for n,wave in enumerate(list(amplitudes.keys())):
        print("INITIALIZE %s"%wave)
        start=False
        if n==0:
            start=True
            WRITE=False
        else:
            WRITE=True
        tsb.rot_write(pert_shape='fft',base_size=nar([size]*3),
                      pert=amplitudes[wave],directory=out_dir_b,
                      k_rot=kint,
                      wave=wave, start=start,write=WRITE, blah=blah)

    tsf = p49_eigen.waves(hx=old_b[0],hy=old_b[1],hz=old_b[2],p=0.6,
                           form='rb96', HydroMethod=4)
    tsf.rot_write(pert_shape='fft',base_size=nar([size]*3),
                  pert=amplitudes['f-'],directory=out_dir_f,
                  k_rot=kint,
                  wave='f-', start=True,write=WRITE, blah=blah)
    tss = p49_eigen.waves(hx=old_b[0],hy=old_b[1],hz=old_b[2],p=0.6,
                           form='rb96', HydroMethod=4)
    tss.rot_write(pert_shape='fft',base_size=nar([size]*3),
                  pert=amplitudes['s+'],directory=out_dir_s,
                  k_rot=kint,
                  wave='s+', start=True,write=WRITE, blah=blah)

if 1:
    field='d'
    print("== fb f- ==")
    print_nz_2(tsb.hat_box['f-'][field])
    print("== fb s+ ==")
    print_nz_2(tsb.hat_box['s+'][field])
    print("== ff f- ==")
    print_nz_2(tsf.hat_box['f-'][field])
    print("== fs s+ ==")
    print_nz_2(tss.hat_box['s+'][field])
if 1:
    sb = p49_plot_tools.chomp(out_dir_b)
    sf = p49_plot_tools.chomp(out_dir_f)
    ss = p49_plot_tools.chomp(out_dir_s)
    out_b="./AAA/test_b_"
    out_s="./AAA/test_s_"
    out_f="./AAA/test_f_"
    analysis={'print_wave':True,
              'plot_fields':0,
              'k_mag':0,
              'k_proj':0}
    print("===== both =====")
    p49_plot_tools.do_stuff(stuff=sb,outdir=out_b,**analysis)
    print("===== fast =====")
    p49_plot_tools.do_stuff(stuff=sf,outdir=out_s,**analysis)
    print("===== slow =====")
    p49_plot_tools.do_stuff(stuff=ss,outdir=out_f,**analysis)
