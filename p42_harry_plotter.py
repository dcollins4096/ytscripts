#nid27562
if 'ef' not in dir():
    execfile('go')
sim_base_dir = {}
sim_base_dir['aq05'] = 'aq05_hydro_m9_32'
sim_base_dir['bq05'] = 'bq05_hydro_m9_32_grav'
sim_base_dir['cq05'] = 'cq05_hydro_m9_32_grav_nodrive'
sim_base_dir['bq04'] = 'bq04_m9_grav_512_L1_J4'
sim_base_dir['aq15'] = 'aq15_ppm_m9_drive2_noamr_128'
sim_base_dir['aq16'] = 'aq16_ppm_m9_drive0_noamr_128'
sim_base_dir['dq04'] = 'dq04_m9_grav_512_L4_J4'
sim_base_dir['cq04'] = 'cq04_m9_nodrive_512_grav'
sim_base_dir['eq04'] = 'eq04_m9_grav_512_L4_J32'
sim_base_dir['aq04'] = 'aq04_hydro_M9_drive_512'
sim_base_dir['aq05'] = 'aq05_hydro_M9_32'
sim_base_dir['aq06'] = 'aq06_hydro_m9_32'
sim_base_dir['aq07'] = 'aq07_hydro_v0.1_32'
sim_base_dir['aq08'] = 'aq08_hydro_m9_stationary'
sim_base_dir['aq09'] = 'aq09_amr'
sim_base_dir['aq10'] = 'aq10_wtf'
sim_base_dir['b02_256'] = '/scratch1/dcollins/Paper08/B2/256/'
sim_base_dir['b2_128'] = 'scratch1/dcollins/Paper08/B2/128/'
sim_base_dir['aq15'] = 'aq15_ppm_m9_drive2_noamr_128'
sim_base_dir['aq16'] = 'aq16_ppm_m9_drive0_noamr_128'
sim_base_dir['aq17'] = 'aq17_ppm_m9_drive1_noamr_32'
sim_base_dir['aq18'] = 'aq18_energy'
sim_base_dir['aq20'] = 'aq20_64_SRGIO_check'
bd = '/scratch/00369/tg456484/Paper42_NewAK'
bd = '/scratch1/dcollins/Paper42_new_turb/'
weight_fields = {'scaled_div_b':'cell_volume'}
methods = {'abs_divb':'mip'}

if 'framelist' not in dir():
    framelist = [76]
fieldlist = ['density','magnetic_field_strength']
fieldlist = ['density'] #,'magnetic_field_strength']
fieldlist = ['density'] #,'TotalEnergy']
fieldlist = ['vorticity_magnitude']
phase_list = [['kinetic_energy','magnetic_energy']]
phase_list = [['density','magnetic_field_strength']]
def square_velocity_divergence(field,data):
    return 4./3*data['velocity_divergence']**2
yt.add_field( 'square_velocity_divergence',function=square_velocity_divergence,units='1/s**2')
def square_velocity_divergence(field,data):
    return np.sign(data['velocity_divergence'])*4./3*data['velocity_divergence']**2
yt.add_field( 'linear_square_velocity_divergence',function=square_velocity_divergence,units='1/s**2',take_log=False)

def forcing_addition(field,data):
    eta = data.get_field_parameter('random_forcing_norm')
    vx,vy,vz = [data['%s-velocity'%vel].in_units('code_velocity').v for vel in 'xyz']
    dx,dy,dz = [data['DrivingField%s'%vel].v for vel in '123']
    d = data['density'].in_units('code_density').v
    return np.abs(d*(eta*(vx*dx+vy*dy+vz*dz) + 0.5*eta*eta*(dx*dx+dy*dy+dz*dz)))
yt.add_field( 'energy_injection', function = forcing_addition) #, take_log=False)



phase_list = [['kinetic_energy','vorticity_squared'],['kinetic_energy','square_velocity_divergence'],['vorticity_squared','square_velocity_divergence']]
phase_list = [['vorticity_squared','square_velocity_divergence']]
phase_list = [['density','square_velocity_divergence']]
phase_list = [['density','vorticity_squared']]
phase_list = [['density','linear_square_velocity_divergence']]
phase_list = [['density','energy_injection']]
phase_list = ['kinetic_energy']
#phase_list+=['velocity_magnitude']
profile_dict = {}
plt.close('all')
#fieldlist = ['DrivingField1','DrivingField2','DrivingField3']
#fieldlist = ['%s-velocity'%s for s in 'xyz']
accumulation = False
project = 'p42'
<<<<<<< /home1/00369/tg456484/yt3_scripts/p42_harry_plotter.py
do_proj = False
do_phase= False
do_profile=False
do_quan = True
if 'quan_list' not in dir():
    quan_list = []
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

for frame_looper in framelist:
    for nsim,sim  in enumerate(simlist):
        if nsim == 1 and False:
            last_time = ds['InitialTime']
            frame = frames[ np.argmin( np.abs(times - last_time)) ]
        else:
            frame = frame_looper


        name  = '%s/%s/DD%04d/data%04d'%(bd,sim_base_dir[sim],frame,frame)
=======
do_proj = False
do_phase= False
do_profile=False
do_injection = False
do_quan = True
if 'quan_list' not in dir():
    quan_list=[]
if do_injection:
    ef('p42_parse_rfout.py')
    frame_list = n
    time_list = t
    eta_list = array[:,26] #VERIFY THIS 
has_prior=False
delta_e = 0.
prior_time=0
prior_ke = 0
delta_t = 1.
prior_norm = 0
print "\n\n\n\n\n"
def quadratic(a,b,c):
    out = -b + np.sqrt(b*b-4*a*c)
    out /= 2*a
    return out
for sim  in ['aq18']:
    if sim in ['b02_256']:
        this_bd = ""
        dirname = 'RS'
        setname = 'restart'
    elif sim in ['b2_128']:
        this_bd = ''
        dirname = 'DD'
        setname = 'data'
    else:
        dirname = 'DD'
        setname = 'data'
        this_bd = bd

    for nf,frame in enumerate(framelist):
        name  = '%s/%s/%s%04d/%s%04d'%(this_bd,sim_base_dir[sim],dirname,frame,setname,frame)
        print name
>>>>>>> /tmp/p42_harry_plotter.py~other.mkKaed
        ds = yt.load(name)
        print name, "t=",ds['InitialTime']
        #print fieldlist
        #continue
        ad=ds.all_data()
<<<<<<< /home1/00369/tg456484/yt3_scripts/p42_harry_plotter.py
        if do_quan:
            quan_list.append( [(ad['density']*ad['%s-velocity'%dim]**2*ad['cell_volume']).sum() for dim in 'xyz'] )
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
=======
        if do_quan:
            #quan_list.append( ad.quantities['WeightedAverageQuantity']('kinetic_energy','cell_volume'))
            #quan_list.append( ad.quantities['Extrema']('density')[1])
            quan_list.append([ds['InitialCycleNumber'],ds['InitialTime']]+ [(ad['density']*ad['%s-velocity'%dim]**2*ad['cell_volume']).sum() for dim in 'xyz'] )

        if do_injection:
            ad.set_field_parameter('random_forcing_norm', np.interp(ds['InitialTime'], time_list,eta_list))
        if do_profile:
            for x_field in  np.unique(nar(phase_list).flatten()):
                z_field = 'cell_mass'
                this_fig, this_ax = profile_dict.get(x_field, plt.subplots(1))
                if not profile_dict.has_key(x_field):
                    profile_dict[x_field] = (this_fig,this_ax)
                #x_field = 'density'
                lims = {'Cooling_Time':[1e12,1e21], 'density':[1e-31,1e-22], 'Temperature':[1e3,1e7]}
                lims['magnetic_field_strength'] = [5e-19,1e-4]
                lims['angular_momentum_magnitude'] = [1e59, 1e71]
                lims['specific_angular_momentum_magnitude'] = [1e25,1e33]
                lims['kinetic_energy'] = [1e-3,1e4]
                lims['velocity_magnitude'] = [0.1,100]
                #if not lims.has_key(y_field):
                #    lims[y_field] = ad.quantities['Extrema'](y_field)
                if not lims.has_key(x_field):
                    lims[x_field] = ad.quantities['Extrema'](x_field)
                phase = yt.create_profile(ad,x_field, #bin_fields=[x_field,y_field], 
                                          fields=[z_field],weight_field=None,
                                          extrema={x_field:lims[x_field]})
                                          #logs={x_field:'False'})
                                          #n_bins=[128,128])
                                          #logs={x_field:x_log,y_field:y_log}) #, n_bins=[32,32])
                labial = '%s_%04d'%(sim,frame)

                color_dict={}
                line_dict={}
                this_ax.plot( 0.5*(phase.x_bins[1:] + phase.x_bins[:-1]), phase['cell_mass'] ,label=labial,c=color_dict.get(sim,'k'),
                             linestyle=line_dict.get(frame,'-'))
                this_ax.set_xlabel(x_field)
                this_ax.set_ylabel(z_field)
                this_ax.set_xscale('log')
                this_ax.set_yscale('log')
                #this_ax.legend(loc=0)
                #plt.clf()
                #plt.plot( 0.5*(phase.x_bins[1:] + phase.x_bins[:-1]), phase['cell_mass'] ,label=labial)
                #plt.savefig('p33y_%s_n%04d'%(sim,frame))

                outname = "%s_%s_%04d_profile_%s_%s.pdf"%(project,sim,frame,x_field,z_field)
                print this_fig.savefig(outname)
                print outname

        if do_phase:
            for x_field, y_field in phase_list:
                #x_field = 'density'
                lims = {'Cooling_Time':[1e12,1e21], 'density':[1e-2,1e2],  'Temperature':[1e3,1e7], 'square_velocity_divergence':[1e-3,1e7],
                        'vorticity_squared':[1e-3,1e7], 'linear_square_velocity_divergence':[-1e5,1e5]}

                if not lims.has_key(y_field):
                    lims[y_field] = ad.quantities['Extrema'](y_field)
                if not lims.has_key(x_field):
                    lims[x_field] = ad.quantities['Extrema'](y_field)
                phase = yt.create_profile(ad,bin_fields=[x_field,y_field],
                                          fields=['cell_mass'],weight_field=None,
                                          extrema={x_field:lims[x_field], y_field:lims[y_field]},
                                          n_bins=[128,128])
                                          #logs={x_field:x_log,y_field:y_log}) #, n_bins=[32,32])
>>>>>>> /tmp/p42_harry_plotter.py~other.mkKaed
                pp = yt.PhasePlot.from_profile(phase)
<<<<<<< /home1/00369/tg456484/yt3_scripts/p42_harry_plotter.py
                pp.set_xlabel(f1)
                pp.set_ylabel(f2)
                print pp.save('%s_%s_n%04d.pdf'%(project,sim,frame))
        if do_proj:
            for ax in 'z':
=======
                pp.set_xlabel(x_field)
                pp.set_ylabel(y_field)
                print pp.save("%s_%s_%04d"%(project,sim,frame))
                if x_field in ['vorticity_squared']:
                    pp.plots[('gas', 'cell_mass')].axes.plot([1e-3,1e7],[1e-3,1e7])
                    print pp.save("%s_%s_%04d"%(project,sim,frame))
        if do_proj:
            for ax in 'x': #yz':
>>>>>>> /tmp/p42_harry_plotter.py~other.mkKaed
                for field in fieldlist:
                    proj = yt.ProjectionPlot(ds,ax,field) #,weight_field=weight_fields.get(field,None),method=methods.get(field,'integrate'))
                    #proj.set_width(50,'Mpc')
                    #proj.annotate_streamlines('Bx','By')
                    #proj.annotate_magnetic_field()
<<<<<<< /home1/00369/tg456484/yt3_scripts/p42_harry_plotter.py
                    if 1:
                        tff = (3*np.pi/(32.*ds['GravitationalConstant']*1))**0.5
                        t_tff = ds['InitialTime']/tff
                        proj.annotate_text([0.05,0.05,0.05],r"$t=%0.2f t_{\rm{ff}}$"%t_tff)
                    if 1:
                        proj.annotate_grids()
=======
                    proj.annotate_grids()
>>>>>>> /tmp/p42_harry_plotter.py~other.mkKaed
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

        if 0:
            vx = ad['x-velocity'].in_units('code_velocity')
            vy = ad['y-velocity'].in_units('code_velocity')
            vz = ad['z-velocity'].in_units('code_velocity')
            d  = ad['density'].in_units('code_density')
            print "Px", np.mean(d*vx)
            print "Py", np.mean(d*vy)
            print "Pz", np.mean(d*vz)
        if 0:
            dv = 1./32**3
            this_time = ds['InitialTime']
            #print prior_ke
            ke = ad['kinetic_energy'].sum().v*dv
            delta_e = ke - prior_ke 
            prior_ke = ke
            delta_t = this_time-prior_time
            prior_time = this_time


            vx = ad['x-velocity'].in_units('code_velocity')
            vy = ad['y-velocity'].in_units('code_velocity')
            vz = ad['z-velocity'].in_units('code_velocity')
            d  = ad['density'].in_units('code_density')
            if nf == 0:
                dvx = ad['x-velocity'].in_units('code_velocity')
                dvy = ad['y-velocity'].in_units('code_velocity')
                dvz = ad['z-velocity'].in_units('code_velocity')
                delta_t = 2.839619e-03 
            if nf > 0:
                stat(vx- vx_pred,'pred vx')
                stat(vy- vy_pred,'pred vy')
                stat(vz- vz_pred,'pred vz')
            stat(vx.v- 0.1, 'm01 vx')
            stat(vy.v- 0.1, 'm01 vy')
            stat(vz.v- 0.1, 'm01 vz')
            print vx


            a = np.sum(0.5*(vx*vx+vy*vy+vz*vz)*d).v
            b = np.sum((vx*dvx+vy*dvy+vz*dvz)*d).v
            a = np.sum(0.5*(dvx*dvx+dvy*dvy+dvz*dvz)*d).v
            c = -ds['RandomForcingEdot'] * delta_t/dv
            norm2=quadratic(a,b,c)
            Edot2 = norm2*b+norm2**2*a+c  #is zero.
            DeltaE2 = np.sum((0.5*((vx+norm2*dvx)*(vx+norm2*dvx)+(vy+norm2*dvy)*(vy+norm2*dvy)+(vz+norm2*dvz)*(vz+norm2*dvz))*d -    0.5*(vx*vx+vy*vy+vz*vz)*d ))*dv
            Edot2 = DeltaE2/delta_t
            vx_pred = vx+ norm2*dvx
            vy_pred = vy+ norm2*dvy
            vz_pred = vy+ norm2*dvz
            print "a %0.2e b %0.2e c %0.2e norm %0.2e Edot2 %0.2d"%(a,b,c,norm2, Edot2)
            print "stuff: cycle %d, dt %0.2e ke %0.5e, prior %0.5e, Edot %0.5e"%( ds['InitialCycleNumber'], delta_t, ke, prior_ke, delta_e/delta_t)



