if 'ef' not in dir():
    execfile('go')
if 'framelist' not in dir():
    framelist = range(2) #range(167) #range(68,112)
    framelist = [92] #range(92)

if 'fieldlist' not in dir():
    fieldlist = []
    fieldlist += ['density']
    #fieldlist += ['magnetic_energy']
#fieldlist += ['DivB','abs_divb'  ]
if 'width' not in dir():
    width=False
ef('p33_sims.py')

def nparticles(this_ds):
    n=0
    for g in this_ds.index.grids:
        n+=sum(g.NumberOfActiveParticles.values())
    return n
if 'simlist' not in dir():
    simlist = ['i02']
for sim in simlist:
    for frame in framelist:
        name  = '%s/%s/DD%04d/DD%04d'%(bd,sim_base_dir[sim],frame,frame)
        print name
        ds = yt.load(name)
        nap = nparticles(ds)
        print "num active partyholes", nap
        ad = ds.all_data()
        for ax in 'x':
            for field in fieldlist:
                proj = yt.ProjectionPlot(ds,ax,field, method=methods.get(field,'integrate'),
                                         weight_field=weight_fields.get(field,None))

                if width:
                    proj.set_width((0.2,'Mpc'))
                #proj.set_zlim('density',1e-8,1e-4)
                #proj.annotate_grids()
                #proj.annotate_streamlines('By','Bz')
                proj.annotate_magnetic_field( scale = 10 )

                proj.set_cmap('density','gray')
                print "NAP", nap
                if nap > 0:
                    proj.annotate_particles(1.0, ptype = 'CenOstriker')
                #stat( ad[field], field)
                outname='p33_%s_n%04d'%(sim,frame)
                print proj.save(outname)
