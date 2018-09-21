if 'ef' not in dir():
    execfile('go')
if 'framelist' not in dir():
    framelist = range(2) #range(167) #range(68,112)
    framelist = [92] #range(92)

if 'fieldlist' not in dir():
    fieldlist = []
    fieldlist += ['density']
    fieldlist = []
    fieldlist += ['magnetic_energy']
#fieldlist += ['DivB','abs_divb'  ]
if 'width' not in dir():
    width=False
ef('p33_sims.py')
if 'do_phase' not in dir():
    do_phase = False
if 'do_proj' not in dir():
    do_proj = True
if 'do_profile' not in dir():
    do_profile = False
def nparticles(this_ds):
    n=0
    for g in this_ds.index.grids:
        n+=sum(g.NumberOfActiveParticles.values())
    return n
if 'simlist' not in dir():
    simlist = ['i02']
if 'phase_list' not in dir():
    #all_fields = all_fields_1
    all_fields = all_fields_simple
    L=len(all_fields)
    phase_list = zip(['density']*L,all_fields)
    phase_list += [['density','magnetic_field_strength']]
    phase_list += [['Temperature','magnetic_field_strength']]
    #phase_list += [['density','Temperature']]

phase_list += [['density','magnetic_field_strength']]
phase_list += [['Temperature','magnetic_field_strength']]
#phase_list = ['magnetic_field_strength'] #,'density','Temperature']
#phase_list = ['density','Temperature']
#phase_list = ['specific_angular_momentum_magnitude']
fr={'ai17':[60, 0]}
fr['ai22'] = fr['ai17']
fr['ai17b'] = fr['ai17']
fr['ai14'] = nar(fr['ai17'])*2
fr['ai16'] = fr['ai14']
fr['ai21'] = fr['ai14']
fr['ai23'] = fr['ai14']
fr['ai24'] = fr['ai14']
color_dict = {'ai16':'r', 'ai17':'g', 'ai21':'b', 'ai22':'k', 'ai23':'r'}
line_dict = {0:'-', 75:'--',150:'--', 5:':', 10:':'}
profile_dict = {}
plt.close('all')
rmap = rainbow_map(len(simlist))
if 'lims' not in dir():
    lims = {'Cooling_Time':[1e12,1e21], 'density':[1e-31,1e-22], 'Temperature':[1e3,1e7]}
    lims['magnetic_field_strength'] = [5e-19,1e-4]
    lims['angular_momentum_magnitude'] = [1e59, 1e71]
    lims['specific_angular_momentum_magnitude'] = [1e25,1e33]
extrema = lims
for nsim,sim in enumerate(simlist):
    #framelist = fr[sim]
    for frame in framelist:
        name  = '%s/%s/DD%04d/DD%04d'%(bd,sim_base_dir[sim],frame,frame)
        ds = yt.load(name)
        print name
        print "TIME %f Myr", ds['InitialTime']
        #nap = nparticles(ds)
        #nap = 0
        #print "num active partyholes", nap
        #continue
        ad = ds.all_data()

        if do_extrema:
            for x_field in  np.unique(nar(phase_list).flatten()):
                this_extrema = ad.quantities['Extrema'](x_field)
                print "%s %s "%(sim,x_field), this_extrema
                if not extrema.has_key(x_field):
                    lims[x_field] = [this_extrema[0].v, this_extrema[1].v]
                else:
                    lims[x_field][0] = min([this_extrema[0].v, lims[x_field][0]])
                    lims[x_field][1] = max([this_extrema[1].v, lims[x_field][1]])
            continue
        ad.set_field_parameter('center',centers.get(sim,nar([0.5,0.5,0.5])))
        #print ad.quantities['Extrema']('specific_angular_momentum_magnitude')
        #continue
        if do_proj:
            for ax in 'xz':
                for field in fieldlist:

                    #proj = yt.ProjectionPlot(ds,ax,field, method=methods.get(field,'integrate'),
                    #j                        weight_field='density') #weight_fields.get(field,None))
                    #L =nar([(10,'kpc'),(10,'kpc'),(25,'kpc')])
                    #R =nar([(50,'kpc'),(50,'kpc'),(35,'kpc')])
                    #C =nar([(30,'kpc'),(30,'kpc'),(30,'kpc')])
                    L = nar([10./60,10./60,25./60])
                    R = nar([50./60,50./60,35./60])
                    C=0.5*(L+R)
                    reg = ds.region(C,L,R)
                    weight_field = weight_fields.get(field,None)
                    p1 = ds.proj(field,ax,data_source=reg,center=C, weight_field=weight_field)
                    proj = p1.to_pw(center=C)

                    if width:
                        proj.set_width((40,'kpc'))
                    #proj.set_zlim('density',1e-8,1e-4)
                    #proj.annotate_grids()
                    #proj.annotate_streamlines('Bx','By')
                    #proj.annotate_magnetic_field()

                    proj.set_cmap('density','gray')
                    #print "NAP", nap
                    #if nap > 0:
                    #    proj.annotate_particles(1.0, ptype = 'CenOstriker')
                    #stat( ad[field], field)
                    outname='p33_%s_n%04d_reg_na'%(sim,frame)
                    print proj.save(outname)
        if do_phase:
            for x_field, y_field in phase_list:
                #x_field = 'density'
                if not lims.has_key(y_field):
                    lims[y_field] = ad.quantities['Extrema'](y_field)
                if not lims.has_key(x_field):
                    lims[x_field] = ad.quantities['Extrema'](y_field)
                phase = yt.create_profile(ad,bin_fields=[x_field,y_field], 
                                          fields=['cell_mass'],weight_field=None,
                                          extrema={x_field:lims[x_field], y_field:lims[y_field]},
                                          n_bins=[128,128])
                                          #logs={x_field:x_log,y_field:y_log}) #, n_bins=[32,32])
                pp = yt.PhasePlot.from_profile(phase)
                pp.set_xlabel(x_field)
                pp.set_ylabel(y_field)
                pp.set_xlim( lims[x_field][0], lims[x_field][1])
                pp.set_ylim( lims[y_field][0], lims[y_field][1])
                print pp.save("%s_%s_%04d"%('p33',sim,frame))
        if do_profile:
            for x_field in  np.unique(nar(phase_list).flatten()):
                z_field = 'cell_mass'
                this_fig, this_ax = profile_dict.get(x_field, plt.subplots(1))
                if not profile_dict.has_key(x_field):
                    profile_dict[x_field] = (this_fig,this_ax)
                #x_field = 'density'
                #lims = {'Cooling_Time':[1e12,1e21], 'density':[1e-31,1e-22], 'Temperature':[1e3,1e7]}
                #if not lims.has_key(y_field):
                #    lims[y_field] = ad.quantities['Extrema'](y_field)
                if not lims.has_key(x_field):
                    lims[x_field] = ad.quantities['Extrema'](x_field)
                phase = yt.create_profile(ad,x_field, #bin_fields=[x_field,y_field], 
                                          fields=[z_field],weight_field=None,
                                          extrema={x_field:lims[x_field]} ) #, y_field:lims[y_field]})
                                          #logs={x_field:'False'})
                                          #n_bins=[128,128])
                                          #logs={x_field:x_log,y_field:y_log}) #, n_bins=[32,32])
                labial = '%s_%04d'%(sim,frame)
                this_ax.plot( 0.5*(phase.x_bins[1:] + phase.x_bins[:-1]), phase['cell_mass'] ,label=labial,c=rmap(nsim),
                             linestyle=line_dict.get(frame,'-'))
                this_ax.set_xlabel(x_field)
                this_ax.set_ylabel(z_field)
                this_ax.set_xscale('log')
                this_ax.set_yscale('log')
                this_ax.legend(loc=0)
                #plt.clf()
                #plt.plot( 0.5*(phase.x_bins[1:] + phase.x_bins[:-1]), phase['cell_mass'] ,label=labial)
                #plt.savefig('p33y_%s_n%04d'%(sim,frame))

                outname = "%s_%s_%04d_profile_%s_%s.pdf"%('p33x',sim,frame,x_field,z_field)
                print this_fig.savefig(outname)
                print outname

if do_profile:
    for x_field in profile_dict.keys():
        this_fig,this_ax = profile_dict[x_field]
        simstr = "%s"*len(simlist)%tuple(simlist)
        outname = "%s_%s_profile_%s_%s.pdf"%('p33_mid',simstr,x_field,z_field)
        this_ax.set_xlabel(x_field)
        this_ax.set_ylabel(z_field)
        this_ax.set_xlim(lims[x_field][0], lims[x_field][1])
        this_ax.set_xscale('log')
        this_ax.set_yscale('log')
        this_ax.legend(loc=0)
        print this_fig.savefig(outname)
        print outname
