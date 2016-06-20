#nid27562
if 'ef' not in dir():
    execfile('go')
framelist = range(5) #range(167) #range(68,112)
sim_base_dir = {}
sim_base_dir['z03'] = 'z03_halo_005036_gas+dm-L4_dedner'
sim_base_dir['e04'] = 'e04_small_verylarge_mhd'
bd = '/mnt/c/scratch/sciteam/dcollins/Paper35_cosmology'
weight_fields = {'scaled_div_b':'cell_volume'}
methods = {'abs_divb':'mip'}

framelist = [10]
fieldlist = ['density','magnetic_field_strength']
fieldlist = ['density'] #,'magnetic_field_strength']
phase_list = [['kinetic_energy','magnetic_energy']]
phase_list = [['density','magnetic_field_strength']]
for sim  in ['e04']:
    for frame in framelist:
        for ax in 'z':
            name  = '%s/%s/RD%04d/RD%04d'%(bd,sim_base_dir[sim],frame,frame)
            continue
            ds = yt.load(name)
            if 0:
                for f1, f2 in phase_list:
                    phase = yt.create_profile(ds.all_data(),bin_fields=[f1,f2], fields=['cell_mass'],weight_field=None)
                    pp = yt.PhasePlot.from_profile(phase)
                    pp.set_xlabel(f1)
                    pp.set_ylabel(f2)
                    print pp.save('p35_%s_n%04d.pdf'%(sim,frame))

            for field in fieldlist:
                proj = yt.ProjectionPlot(ds,ax,field) #,weight_field=weight_fields.get(field,None),method=methods.get(field,'integrate'))
                proj.set_width(50,'Mpc')
                #proj.annotate_streamlines('Bx','By')
                proj.annotate_magnetic_field()
                print proj.save('p35_%s_n%04d_field.pdf'%(sim,frame))

