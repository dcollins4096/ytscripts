#nid27562
if 'ef' not in dir():
    execfile('go')
framelist = range(5) #range(167) #range(68,112)
sim_base_dir = {}

bd = '/scratch/00369/tg456484/Paper42_NewAK'
bd = '/scratch1/dcollins/Paper42_new_turb/'
bd = '/home/dcollins4096/ASTRO_12_3/runs'
weight_fields = {'scaled_div_b':'cell_volume'}
methods = {'abs_divb':'mip'}

framelist = range(0,240,10) + [236] 
fieldlist = ['density','magnetic_field_strength']
fieldlist = ['density'] #,'magnetic_field_strength']
phase_list = [['kinetic_energy','magnetic_energy']]
phase_list = [['density','magnetic_field_strength']]
#fieldlist = ['DrivingField1','DrivingField2','DrivingField3']
framelist = range(0,500,100) + [418]
#fieldlist = ['%s-velocity'%s for s in 'xyz']
accumulation = False
project = 'p10'
do_proj = True
do_phase= False
do_profile=False
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
field_list = ['density']
framelist = [0,10,600]

for sim  in [190]:
    if sim not in sim_base_dir:
        sim_base_dir[sim]=sim

    dirname = 'DD'
    setname = 'data'
    for nf,frame in enumerate(framelist):
        name  = '%s/%s/%s%04d/%s%04d'%(bd,sim_base_dir[sim],dirname,frame,setname,frame)
        print name
        ds = yt.load(name)
        ad=ds.all_data()
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
            for ax in 'x': #yz':
                for field in fieldlist:
                    proj = yt.ProjectionPlot(ds,ax,field) #,weight_field=weight_fields.get(field,None),method=methods.get(field,'integrate'))
                    #proj.set_width(50,'Mpc')
                    #proj.annotate_streamlines('Bx','By')
                    #proj.annotate_magnetic_field()
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



