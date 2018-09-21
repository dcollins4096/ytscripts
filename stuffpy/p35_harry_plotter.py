#nid27562
import scatter_fit
if 'ef' not in dir():
    execfile('go')
framelist = range(5) #range(167) #range(68,112)
sim_base_dir = {}
sim_base_dir['z03'] = 'z03_halo_005036_gas+dm-L4_dedner'
sim_base_dir['e04'] = 'e04_small_verylarge_mhd'
bd = '/mnt/c/scratch/sciteam/dcollins/Paper35_cosmology'
bd = '/Users/dcollins/scratch/P33P35/'
sim_base_dir['h005016'] = ''
bd1 = '/mnt/c/scratch/sciteam/dcollins/Paper35_cosmology'
bd2 = '/u/sciteam/brittons/projects/bw_mw/mean/halo_005016/gas+dm-L4/'
bd = bd1
weight_fields = {'scaled_div_b':'cell_volume'}
methods = {'abs_divb':'mip'}

framelist = [1]
fieldlist = ['density','magnetic_field_strength']
#fieldlist = [] #
phase_list = [['kinetic_energy','magnetic_energy']]
phase_list = [['density','magnetic_field_strength']]
for sim  in ['z03']:
    for frame in framelist:
        for ax in 'z':
            setprefix = 'DD'
            name  = '%s/%s/%s%04d/%s%04d'%(bd,sim_base_dir[sim],setprefix,frame,setprefix,frame)
            ds = yt.load(name)
            ad=ds.all_data()
            if 1:
                for f1, f2 in phase_list:
                    phase = yt.create_profile(ad,bin_fields=[f1,f2], fields=['cell_mass'],weight_field=None)
                    profile = yt.create_profile(ad, f1, f2, weight_field='cell_mass')
                    pp = yt.PhasePlot.from_profile(phase)
                    pp.set_xlabel(f1)
                    pp.set_ylabel(f2)
                    print pp.save('p35_%s_n%04d.pdf'%(sim,frame))
                    this_axes = pp.plots[('gas', 'cell_mass')].axes
                    #pp.plots[('gas', 'cell_mass')].axes.plot([1e-32,1e-26],[1e-12,1e-6])
                    xb = 0.5*(profile.x_bins[1:]+profile.x_bins[0:-1])
                    field_strength = profile[f2]
                    if 1:
                        plt.clf()
                        plt.plot(xb,field_strength,marker='x')
                        plt.xscale('log'); plt.yscale('log')
                    fit = scatter_fit.scatter_fit(this_axes,xb,field_strength) #, fit_range=[1e-31,1e-28])
                    powerline(this_axes, fit['x'][0], fit['x'][1], fit['y'][0],0.66,{'c':'y'})
                    plt.savefig('herp.png')
                    this_axes.plot(xb, field_strength, marker='x')
                    print pp.save('p35_%s_n%04d'%(sim,frame))

            for field in fieldlist:
                proj = yt.ProjectionPlot(ds,ax,field) #,weight_field=weight_fields.get(field,None),method=methods.get(field,'integrate'))
                #proj.set_width(50,'Mpc')
                #proj.annotate_streamlines('Bx','By')
                proj.annotate_magnetic_field()
                print proj.save('p35_%s_n%04d_field'%(sim,frame))

