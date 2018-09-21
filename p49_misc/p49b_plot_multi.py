
from go_lite_pyt3 import *
import yt
from yt.funcs import mylog
mylog.setLevel(50)
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
from p49_print_tools import *
plt.close('all')
#frame_list=range(0,60,10)
frame_list=[0,50]
this_formt = 'png'
get_from='yt'
plt.close('all')
wave_list  =['f-', 'a-','s-','c','f+','a+','s+']
if 'readit' not in dir():
    readit=True
if readit:
    p1d={}
    if 'lazy_ds' not in dir() or True:
        lazy_ds = {}
    runs = ['rx11_mixed_s+', 'rx11_mixed_1','rx11_mixed_f-']
    od = "%s_%s_%s"%tuple(runs)
    data_dir_base = '/Users/dcollins/scratch/Paper49b_play/Spectral/'
    outdir = '/Users/dcollins/RESEARCH2/Paper49_EBQU/2018-06-12-p49b/EigenTests_PostFFT/'+od
    print(outdir)
    data_dirs = ["%s/%s"%(data_dir_base,r) for r in runs]
    for frame in frame_list: 
        p1d[frame]={}
        for nr,r in enumerate(runs):
            ds_fname = "%s/DD%04d/data%04d"%(data_dirs[nr],frame,frame)
            print("reading %s"%ds_fname)
            ds = yt.load(ds_fname)
            stuff = p49_eigen.get_cubes_cg(ds)
            p1d[frame][r]=p49_plot_tools.plot_wave_mag(stuff=stuff, output_name = "AAA/tmp.png")

if 1:
    for frame in frame_list: 
        Nx = 1 #len(runs)
        Xinch = Nx*4
        Ny = len(wave_list)
        Yinch = Ny*4

        fig = plt.figure(figsize=(Xinch,Yinch)) 
        ax = [fig.add_subplot(Ny,Nx,i+1) for i in range(Ny*Nx)]

        for a in ax:
            a.set_xticklabels([])
            a.set_yticklabels([])
            a.set_aspect('equal')

        fig.subplots_adjust(wspace=0, hspace=0)

        oots=[None]*len(ax)

        count=-1
        for nw,wave in enumerate(wave_list):
            count+=1
            for r in runs:
                ax[count].plot(p1d[frame][r]['bin_center'],p1d[frame][r][wave],label=r)
                ax[count].set_yscale('log')
                ax[count].set_ylim([1e-10,1e-7])
            ax[count].legend(loc=0)
        oname = "%s/multi_amp_%04d"%(outdir,frame)
        fig.savefig(oname)


            



        #,stuff=None,output_name="thig.png", nbins=None)

