
from go_lite_pyt3 import *
import yt
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
from p49_stuff import *
plt.close('all')
frame_list=[0]#,50]#range(0,60,10)
this_formt = 'png'
get_from='yt'
if 1:
    if 'lazy_ds' not in dir() or True:
        lazy_ds = {}
    for frame in frame_list: 
        if 0:
            this_name = 'r801'
            directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/test_new'
            #directory = '/Users/dcollins/scratch/Paper49b_play/Spectral/test_good'
        if 0:
            outdir = '/Users/dcollins/RESEARCH2/Paper49_EBQU/2018-06-12-p49b/EigenTests/r701'
            directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/701_repeat'
            prefix = 'r701_orig'
        if 0:
            outdir = '/Users/dcollins/RESEARCH2/Paper49_EBQU/2018-06-12-p49b/EigenTests_PostFFT'
            run='r701_new'
        directory_e = '/Users/dcollins/scratch/Paper49b_play/Eigen/'
        directory_s = '/Users/dcollins/scratch/Paper49b_play/Spectral/'
        directory=directory_e
        #run = 'r701_new'
        #run = 'r701_spec'; directory=directory_s
        #run = 'r701b_hydro4'; 
        #run = 'r701c_32'; directory=directory_s
        #run = 'rx01_202'; directory=directory_s
        #run = 'rx02_-202'; directory=directory_s
        #run = 'rx03_220'; directory=directory_s
        #run = 'rx04_221'; directory=directory_s
        #run = 'rx05_212'; directory=directory_s
        #run = 'rx06_212_b2'; directory=directory_s
        #run = 'rx07_212_128'; directory=directory_s
        run = 'rx07b_212_128_HLLD'; directory=directory_s

        directory = "%s/%s"%(directory,run)
        outdir = '/Users/dcollins/RESEARCH2/Paper49_EBQU/2018-06-12-p49b/EigenTests_PostFFT/'+run
        #outdir = "./AAA/%s"%run
        outdir += "_repeat"
        prefix=run
        if True: #frame not in lazy_ds
            if get_from=='yt':
                ds_fname = "%s/DD%04d/data%04d"%(directory,frame,frame)
                print("reading %s"%ds_fname)
                ds = lazy_ds.get(frame,yt.load(ds_fname))
                stuff = p49_eigen.get_cubes_cg(ds)
                stuff['ds']=ds
                lazy_ds[frame]=stuff

        def maxis(inarr,dim):
            arr=np.zeros_like(inarr)
            ok =  inarr>1e-13
            arr[ok] = inarr[ok]
            out = np.log(np.max(np.abs(arr),axis=dim))

            #out = np.sum(np.abs(inarr),axis=dim)
            #out[ out==0] == -1
            return out
        analysis={'print_wave':False,
                  'plot_fields':1,
                  'k_mag':1,
                  'k_proj':1,
                  'k_func':maxis
                 }
        print("Dumping into %s/%s_%04d"%(outdir,prefix,frame))
        p49_plot_tools.do_stuff(stuff=stuff,outdir="%s/%s_%04d_"%(outdir,prefix,frame),**analysis)
        #p1d=p49_plot_tools.plot_wave_mag(stuff=stuff,output_name="%s/%s_%04d"%(outdir,prefix,frame))

