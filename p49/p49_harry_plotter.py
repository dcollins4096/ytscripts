#nid27562
if 'ef' not in dir():
    execfile('go')
sim_base_dir = {}

sim_base_dir['aw05'] = 'aw05_M1_MA1_1024'
sim_base_dir['aw06'] = 'aw06_M1_MA0.3_1024'

bd = '/scratch/00369/tg456484/Paper49_EB/'
weight_fields = {'scaled_div_b':'cell_volume'}
methods = {'abs_divb':'mip'}

if 'framelist' not in dir():
    framelist = [0,10,20,30,40]
fieldlist = ['density','magnetic_field_strength']
fieldlist = ['density'] #,'magnetic_field_strength']
fieldlist = ['density'] #,'TotalEnergy']
phase_list = [['kinetic_energy','magnetic_energy']]
phase_list = [['density','magnetic_field_strength']]
accumulation = False
project = 'p49'
do_proj = True
do_phase= False
do_profile=False

if 0:
    try:
        fff = open('/scratch/00369/tg456484/Paper42_NewAK/aq16_ppm_m9_drive0_noamr_128/OutputLog')
        names = []
        times = []
        frames = []
        for l in fff.readlines():
            nw= no_whites(l.split(" "))
            times.append(float(nw[4]))
            names.append(nw[2])
            frames.append( int(names[-1][4:8]) )
        times=nar(times)
    except:
        raise

if 'simlist' not in dir():
    simlist = ['aw05','aw06']

for frame_looper in framelist:
    for nsim,sim  in enumerate(simlist):
        if nsim == 1 and False:
            last_time = ds['InitialTime']
            frame = frames[ np.argmin( np.abs(times - last_time)) ]
        else:
            frame = frame_looper


        name  = '%s/%s/DD%04d/data%04d'%(bd,sim_base_dir[sim],frame,frame)
        ds = yt.load(name)
        print name, "t=",ds['InitialTime']
        #print fieldlist
        #continue
        ad=ds.all_data()
        if 0:
            px = ad.quantities['WeightedAverageQuantity']('momentum_x','cell_volume')
            py = ad.quantities['WeightedAverageQuantity']('momentum_y','cell_volume')
            pz = ad.quantities['WeightedAverageQuantity']('momentum_z','cell_volume')
            print px
        if do_profile:
            field = 'density'
            prof = yt.create_profile(ds.all_data(),field,fields='cell_mass',weight_field=None,accumulation=accumulation,
                                    fractional=True)
            plt.clf()
            plt.plot(0.5*(prof.x_bins[1:]+prof.x_bins[0:-1]),prof['cell_mass'], label=sim)
            plt.xscale('log'); plt.yscale('log')
            profname = '%s_prof_%s_n%04d.pdf'%(project, field, frame)
            plt.savefig(profname)
            print profname
        if do_phase:
            for f1, f2 in phase_list:
                phase = yt.create_profile(ds.all_data(),bin_fields=[f1,f2], fields=['cell_mass'],weight_field=None)
                pp = yt.PhasePlot.from_profile(phase)
                pp.set_xlabel(f1)
                pp.set_ylabel(f2)
                print pp.save('%s_%s_n%04d.pdf'%(project,sim,frame))
        if do_proj:
            for ax in 'z':
                for field in fieldlist:
                    proj = yt.ProjectionPlot(ds,ax,field) #,weight_field=weight_fields.get(field,None),method=methods.get(field,'integrate'))
                    #proj.set_width(50,'Mpc')
                    #proj.annotate_streamlines('Bx','By')
                    #proj.annotate_magnetic_field()
                    if 1:
                        tff = (3*np.pi/(32.*ds['GravitationalConstant']*1))**0.5
                        t_tff = ds['InitialTime']/tff
                        proj.annotate_text([0.05,0.05,0.05],r"$t=%0.2f t_{\rm{ff}}$"%t_tff)
                    if 1:
                        proj.annotate_grids()
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



