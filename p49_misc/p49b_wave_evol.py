


if 'ef' not in dir():
    execfile('go')
    for i in range(3):
        print("====================")
import enzo_write
reload(enzo_write)
import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
import matplotlib.colors as colors



def nz(field):
    nz = np.abs(field) > 1e-13
    return field[nz]
frame_list=[0]
this_formt = 'png'
get_from='ic'
#plot_style = 'r_theta'
#plot_style = 'proj'
plot_style = 'hist_tot'


if 1:
    if 'lazy_ds' not in dir():
        lazy_ds = {}
    for frame in frame_list: 
        if 1:
            if 0:
                this_name = 'y701'
                directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/y701_rb96_fft_f-_play'
            if 0:
                this_name = 'r801'
                directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/r801_rj95_110_f-'
            if 0:
                this_name = 'rA01'
                directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/rA01_rb96_110_f-'
            if 1:
                this_name = 'rB01'
                directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/rB01_rb_several'
            #frame//wave//xy//xz//yz?//real//imag//magnitude//phase
            #https://matplotlib.org/users/colormapnorms.html
            plt.close('all')
            if frame not in lazy_ds:
                if get_from=='yt':
                    ds = lazy_ds.get(frame,yt.load("%s/DD%04d/data%04d"%(directory,frame,frame)))
                    stuff = p49_eigen.get_cubes_cg(ds)
                    #lazy_ds[frame]=stuff
                elif get_from=='ic':
                    this_name = 'rB01_ic'
                    stuff = p49_plot.tools.chomp(directory=directory)
                else:
                    print("Extant Stuff")
                #lazy_ds[frame]=ds
            else:
                stuff = lazy_ds[frame]
            print_fields = False
            print_waves = True
            #these_means = stuff['means']
            #these_ffts  = p49_eigen.get_ffts(stuff['cubes'], these_means)
            #kall,wut=p49_eigen.rotate_back(these_ffts, these_means)
            #kmag = (kall[0,...]**2+kall[1,...]**2+kall[2,...]**2)**0.5
        if plot_style == 'hist_tot':
            oname = '%s_%04d_hist.%s'%(this_name, frame, this_formt)
            p49_plot_tools.plot_wave_mag(stuff=stuff,output_name=oname) 



            """
                fig = plt.figure(figsize=(8,8)) # Notice the equal aspect ratio
                fig.suptitle('%s_%04d %s'%(this_name,frame,wave))
                #ax = [fig.add_subplot(1,1,i+1) for i in range(6)]
                ax = [fig.add_subplot(1,1,1,projection='polar')]

                for a in ax:
                    a.set_xticklabels([])
                    a.set_yticklabels([])
                    a.set_aspect('equal')


                all_angle = np.angle(this_fft)
                flag = np.abs(this_fft) > 1e-9
                this_kmag = kmag[flag]
                this_angle = all_angle[flag]
                oname = '%s_%04d_%s_rtheta.%s'%(this_name, frame, wave, this_formt)
                ax[0].scatter(this_angle, this_kmag)
                for a in ax:
                    a.set_rmax(16)
                fig.savefig(oname)
                print(oname)
                """
        if plot_style == 'r_theta':
            p49_plot_tools.plot_k_rad(wut=wut,prefix="%s_%04d"%(this_name,frame))


        if plot_style == 'proj':
            p49_plot_tools.plot_k_proj(wut=wut,prefix="%s_%04d"%(this_name,frame))

if 0:
    #old shit?
    #Test.  Frame 0 has only f-.
    frame = 0
    directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/y701_rb96_fft_f-_play'
    ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
    stuff = p49_eigen.get_cubes_cg(ds)
    these_means = stuff['means']
    these_ffts  = p49_eigen.get_ffts(stuff['cubes'], these_means)
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
            if ( np.abs(thisthing)>1e-12).sum() > 0:
                bang_or_not = "!!!"*8 + " meann %0.2e max %0.2e"%(np.mean(np.abs(thisthing)),np.abs(thisthing).max())
            print("=== Wave %s %s"%(wave, bang_or_not))
            s1 = str(nz(thisthing.real).size)
            s2 = str(nz(thisthing.imag).size)
            print("wave real nz %s imag %s"%(s1,s2))
