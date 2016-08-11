#nid27562
if 'ef' not in dir():
    execfile('go')
framelist = range(5) #range(167) #range(68,112)
sim_base_dir = {}
sim_base_dir['aq05'] = 'aq05_hydro_m9_32'
bd = '/scratch/00369/tg456484/Paper42_NewAK'
weight_fields = {'scaled_div_b':'cell_volume'}
methods = {'abs_divb':'mip'}

framelist = [183]
fieldlist = ['density','magnetic_field_strength']
fieldlist = ['density'] #,'magnetic_field_strength']
phase_list = [['kinetic_energy','magnetic_energy']]
phase_list = [['density','magnetic_field_strength']]
accumulation = False
project = 'p42'
proj = True
phase= False
profile=False
for sim  in ['aq05']:
    for frame in framelist:
        name  = '%s/%s/DD%04d/data%04d'%(bd,sim_base_dir[sim],frame,frame)
        print name
        ds = yt.load(name)
        ad=ds.all_data()
        if profile:
            field = 'density'
            prof = yt.create_profile(ds.all_data(),field,fields='cell_mass',weight_field=None,accumulation=accumulation,
                                    fractional=True)
            plt.clf()
            plt.plot(0.5*(prof.x_bins[1:]+prof.x_bins[0:-1]),prof['cell_mass'], label=sim)
            plt.xscale('log'); plt.yscale('log')
            profname = '%s_prof_%s_n%04d.pdf'%(project, field, frame)
            plt.savefig(profname)
            print profname
        if phase:
            for f1, f2 in phase_list:
                phase = yt.create_profile(ds.all_data(),bin_fields=[f1,f2], fields=['cell_mass'],weight_field=None)
                pp = yt.PhasePlot.from_profile(phase)
                pp.set_xlabel(f1)
                pp.set_ylabel(f2)
                print pp.save('%s_%s_n%04d.pdf'%(project,sim,frame))
        if proj:
            for ax in 'z':
                for field in fieldlist:
                    proj = yt.ProjectionPlot(ds,ax,field) #,weight_field=weight_fields.get(field,None),method=methods.get(field,'integrate'))
                    #proj.set_width(50,'Mpc')
                    #proj.annotate_streamlines('Bx','By')
                    #proj.annotate_magnetic_field()
                    print proj.save('%s_%s_n%04d'%(project,sim,frame))
        if 0:
            vx = ad['x-velocity'].in_units('code_velocity')
            vy = ad['y-velocity'].in_units('code_velocity')
            vz = ad['z-velocity'].in_units('code_velocity')
            d  = ad['density'].in_units('code_density')
            vx_bar = np.mean(vx)
            vy_bar = np.mean(vy)
            vz_bar = np.mean(vz)
            v2 =  (vx-vx_bar)**2 + (vy-vy_bar)**2 + (vz-vz_bar)**2 
            sigma2 = np.mean(v2)
            ke = (d*(v2) * ad['cell_volume'].in_units('code_length**3')).sum()
            if 'ke_sequence' not in dir():
                ke_sequence = {}
            cycle = ds['InitialCycleNumber']
            time = ds['InitialTime']
            ke_sequence[cycle] = [time, ke, sigma2]



