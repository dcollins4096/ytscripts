
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
blah=False
if 0:
    #setup FFT init
    blah=False
    #Target woodshed
    #k3 = p49_eigen.make_k_freqs(size)
    ks = p49_eigen.make_k_freqs_and_int(size)
    k3 = ks['k_freq']
    kint = ks['k_int']
    kint_norm = np.sqrt( kint[0,...]**2 + kint[1,...]**2 + kint[2,...]**2 )
    k_norm = np.sqrt( k3[0,...]**2 + k3[1,...]**2 + k3[2,...]**2 )
    ok = k_norm > 0


amplitude=1e-3
ampl=np.zeros_like(k_norm)*(1+0j)
old_b   = nar([1.0, 1.41421, 0.5])
if 0:
    #single wave.  Works, man.
    ampl[2,0,2]=amplitude
    directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/rx01_202'

if 0:
    #reverse wave.  Works.
    ampl[-2,0,2]=amplitude
    directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/rx02_-202'

if 0:
    #single wave, different angle.  Works fine.
    ampl[2,2,0]=amplitude
    directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/rx03_220'

if 0:
    #single wave, more angle.  Generates even odd.  fascinating.
    #Seems to get corrugated in Vz and Bz first.
    ampl[2,2,1]=amplitude
    directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/rx04_221'

if 0:
    #single wave, more angle, rotated.  Hx gets corrugated first?
    ampl[2,1,2]=amplitude
    directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/rx05_212'

if 1:
    #single wave, more even-odd hunting.  Even-odd in y? Yes, probably.
    old_b   = nar([1.0, 0.5, 1.41421])
    ampl[2,1,2]=amplitude
    directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/rx06_212_b2'

    #kludge to track down a typo somewhere
    #directory= '/Users/dcollins/scratch/Paper49b_play/Spectral/rx08_k-53'

if 0:
    if size != 128:
        print("SIZE ERROR THIS ONE NEEDS TO BE LARGE")
        print("SIZE ERROR THIS ONE NEEDS TO BE LARGE")
        print("SIZE ERROR THIS ONE NEEDS TO BE LARGE")
        print("SIZE ERROR THIS ONE NEEDS TO BE LARGE")
    #single wave, more even-odd hunting.  
    #WITH RiemannSolver=1, this is a horrid even-odd. 
    #RiemannSolver =6 works nicely.
    old_b   = nar([1.0, 0.5, 1.41421])
    ampl[2,1,2]=amplitude
    directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/rx07_212_128'
    directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/rx07b_212_128_HLLD'
else:
    if size != 32:
        print("SIZE ERROR needs to be 32")
        print("SIZE ERROR needs to be 32")
        print("SIZE ERROR needs to be 32")
        print("SIZE ERROR needs to be 32")


if 0:
    #like r701.  Good.
    ampl=np.zeros_like(k_norm)*(1+0j)
    wave='f-'
    ts_temp = p49_eigen.waves(hx=1.0,hy=1.41421,hz=0.5,p=0.6,this_wave=wave, form='rb96')
    ratio = ts_temp.speeds['cf']/ts_temp.speeds['aa']
    ampl[1,0,0]=1e-6*ratio
    directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/r701c_32'


if 0:
    #ensemble In progress
    theta = np.pi*2*np.random.random(ok.sum())
    phase = 1.# np.exp(theta*1j)
    ampl[ok] = k_norm[ok]**(-5./3)*phase
    mask = k3[0,...]>0
    mask = np.logical_and(mask,k3[1,...]>0)
    mask = np.logical_and(mask,k3[2,...]>0)
    bmask = mask == False
    ampl[bmask]=0.
    ampl[:]=0.
    ampl[1,5,7] = amplitude
directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/rx08_k-53'

if 0:
    #test the constituents
    import matplotlib.colors as colors
    img_args={ 'origin':'lower','interpolation':'nearest'}
    img_args['norm'] = colors.Normalize(vmin=0,vmax=np.abs(ampl).max())
    p49_plot_tools.six_plot(ampl,oname = 'AAA/ampl.png',img_args=img_args)
    p49_plot_tools.six_plot(k_norm,oname = 'AAA/k_norm.png')
    p49_plot_tools.six_plot(kint_norm,oname = 'AAA/k_int.png')

if 0:
    ampl=np.ones_like(k_norm)*(amplitude+0j)
    this_mask = np.ones_like(ampl,dtype=bool)
    this_mask[-5,7,1]=True
    blank_mask = this_mask == False
    ampl[blank_mask]=0.


if 0:
    #works.
    directory_good = '/Users/dcollins/scratch/Paper49b_play/Spectral/test_good'
    kint_good = nar([[2,1],[0,0],[2,0]])
    ampl_good = nar([amplitude,0])
    blah_good=True
if 0:
    #works.
    kint_good = nar([[-2,1],[2,0],[1,0]])
    ampl_good = nar([amplitude,0])
    blah_good=True

if 'wave' not in dir():
    wave = 's+'

if 'unlazy' not in dir():
    unlazy=True
if 'ts1' not in dir() or unlazy:
    WRITE=True

    if 'ampl_good' in dir():
        ts1 = p49_eigen.waves(hx=old_b[0],hy=old_b[1],hz=old_b[2],p=0.6,
                              this_wave=wave, form='rb96', HydroMethod=4)
        ts1.rot_write(pert_shape='fft',base_size=nar([size]*3),
                      pert=ampl_good,directory=directory_good,
                      k_rot=kint_good,
                      wave=wave, start=True,write=WRITE, blah=blah_good)

    if 'directory' not in dir():
        directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/test_new'
    ts2 = p49_eigen.waves(hx=old_b[0],hy=old_b[1],hz=old_b[2],p=0.6,
                          this_wave=wave, form='rb96', HydroMethod=4)
    print("wrtie the what?")
    ts2.rot_write(pert_shape='fft',base_size=nar([size]*3),
                  pert=ampl,directory=directory,
                  k_rot=kint,
                  wave=wave, start=True,write=WRITE, blah=blah)

directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/rx08_k-53'
sn = p49_plot_tools.chomp(directory)

analysis={'print_wave':True,
          'plot_fields':0,
          'k_mag':0,
          'k_proj':0}
p49_plot_tools.do_stuff(stuff=sn,outdir="%s/test_new"%od,**analysis)
