
from go_lite_pyt3 import *
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
from p49_print_tools import *

size = 32
if 1:
    ks = p49_eigen.make_k_freqs_and_int(size)
    k3 = ks['k_freq']
    kint = ks['k_int']

if 0:

    tn = ts1
    to = t_works
    c1 = t_works.cubes
    c2 = ts1.cubes
    h1 = t_works.all_hats
    h2 = ts1.all_hats

    def ddd(field):
        a = t_works.cubes[field]-ts1.cubes[field]
        return a

    nza 

if 1:
    b1='/Users/dcollins/scratch/Paper49b_play/Spectral/test_good'
    b2='/Users/dcollins/scratch/Paper49b_play/Spectral/test_new'
    sim= ""
    sim = 'r401_rj95_sq_f-'
    b1 = '/Users/dcollins/scratch/Paper49b_play/Eigen'
    b2 = '/Users/dcollins/scratch/Paper49b_play/Eigen_Take1'
    b1 = '/Users/dcollins/scratch/Paper49b_play/Spectral/';sim1='rx08_k-53'
    b2 = '/Users/dcollins/scratch/Paper49b_play/Spectral/';sim2='rx06_212_b2'
    #sim = 'r501_rj95_fft_f-'
    #sim1 = 'r701_rb96_fft_f-'
    #sim2 = 'r701_rb96_fft_f-'

    file_list = ["Bx_16.h5", "By_16.h5", "Bz_16.h5", "density_16.h5", "x-velocity_16.h5", "y-velocity_16.h5", "z-velocity_16.h5"]
    file_list += ['GasPressure_16.h5']
    #file_list += ['TotalEnergy_16.h5']

    for f in file_list:
        fname1 = "%s/%s/%s"%(b1,sim1,f)
        fp1 = h5py.File(fname1,'r')
        fname2="%s/%s/%s"%(b2,sim2,f)
        fp2 = h5py.File(fname2,'r')
        try:
            a1=fp1[f][:]
            a2=fp2[f][:]
            print("%s %20s %0.2e"%(sim,f,np.abs(a1-a2).sum()/np.sum(a1)))
        except:
            raise
        finally:
            fp1.close()
            fp2.close()

if 0:
    sim = 'r701_rb96_fft_f-'
    s1 = p49_plot_tools.chomp("%s/%s"%(b2,sim),HydroMethod=6)
#p49_plot_tools.do_stuff(stuff=s1,plot_fields=True,outdir="./p49b_r701_orig_")
    s2 = p49_plot_tools.chomp("%s/%s"%(b1,sim),HydroMethod=6)
#p49_plot_tools.do_stuff(stuff=s2,plot_fields=True,outdir="./p49b_r701_broke_")
    for field in s1['cubes'].stuff.keys():
        f1 = s1['cubes'][field]
        nz=f1!=0
        f2 = s2['cubes'][field]
        #dd = np.sum(np.abs(s1['cubes'][field]/s2['cubes'][field]))
        #print("%s %0.2e"%(field, dd  ))
        print("%s 1max %0.2e 2max %0.2e rat %0.2e"%(field,f1.max(),f2.max(),f1.max()/f2.max()))
