
# general numbers:
# n = 0.5 cm^-3; B=6e-6 G; T=7500K gamma=1.666 mu=1.3


if 'ef' not in dir():
    execfile('go')
framelist = range(5) #range(167) #range(68,112)
sim_base_dir = {}
sim_base_dir['a03'] = 'a03_long_subalfv'
bd = '/scratch1/dcollins/Paper42_new_forcing/'
weight_fields = {'scaled_div_b':'cell_volume'}
methods = {'abs_divb':'mip'}

framelist = [100]
fieldlist=['density']
phase_list=[['density','magnetic_field_strength']]
for sim  in ['a03']:
    for frame in framelist:
        for ax in 'z':
            name  = '%s/%s/DD%04d/data%04d'%(bd,sim_base_dir[sim],frame,frame)
            ds = yt.load(name)
            ad=ds.all_data()
            if 1:
                for f1, f2 in phase_list:
                    phase = yt.create_profile(ds.all_data(),bin_fields=[f1,f2], fields=['cell_mass'],weight_field=None)
                    pp = yt.PhasePlot.from_profile(phase)
                    pp.set_xlabel(f1)
                    pp.set_ylabel(f2)
                    print pp.save('p49_%s_n%04d.pdf'%(sim,frame))

            for field in fieldlist:
                proj = yt.ProjectionPlot(ds,ax,field) #,weight_field=weight_fields.get(field,None),method=methods.get(field,'integrate'))
                #proj.set_width(50,'Mpc')
                proj.annotate_streamlines('Bx','By')
                proj.annotate_magnetic_field()
                print proj.save('p49_%s_n%04d_field.pdf'%(sim,frame))

