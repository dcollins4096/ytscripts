if 'ef' not in dir():
    execfile('go')
def tracer_number_density(field,data):
    norm = data.get_field_parameter('particle_norm')

    cic = data[('deposit','all_cic')].in_units('code_density')
    number_density = cic/norm
    return number_density
these_units='1/code_length**3'
#these_units = 'code_density'
yt.add_field('tracer_number_density',function=tracer_number_density,# units=these_units,
            validators=[yt.ValidateParameter("particle_norm")])

if 0:
    basedir=[]
    name=[]
    prof=[]
    labels=[]
    ds=[]
    basedir.append('/scratch1/dcollins/Paper19/u07-64-driving-nograv');name.append('u07');labels.append('Trilinear')
    basedir.append('/scratch1/dcollins/Paper19/u11_u07_newtracer');name.append('u11');labels.append('Flux')
    basedir.append('/scratch1/dcollins/Paper19/u12_montecarlo');name.append('u12');labels.append('MC')
    plt.clf()
    offsets = [1.0,1.01,1.1]
    for n, nm,bd in zip(range(len(name)),name,basedir):
        fname = '%s/DD%04d/data%04d'%(bd,frame,frame)
        ds.append(yt.load(fname))
        ad=ds[-1].all_data()
        ad.set_field_parameter('particle_norm',8*ds[-1].index.grids[0]['particle_mass'][0]/(ds[-1].index.grids[0].dds.prod()))
        if n == 0:
            prof = yt.ProfilePlot(ad,'density','cell_volume',weight_field=None,label=r'$\rho/\rho_0$')
            the_x = prof.profiles[0].x_bins
            xxx = 0.5*(the_x[1:]+the_x[0:-1])
            plt.plot(np.log(xxx),np.log10(prof.profiles[0]['cell_volume']),label=r'$\rho/\rho_0$')
        prof = yt.ProfilePlot(ad,'tracer_number_density','cell_volume',weight_field=None)
        the_x = prof.profiles[0].x_bins
        xxx = 0.5*(the_x[1:]+the_x[0:-1])
        volume = prof.profiles[0]['cell_volume']
        plt.plot(np.log(xxx),np.log10(volume),label=labels[n])
        print (xxx*volume).sum()
    plt.xlabel('log density')
    plt.legend(loc=0)
    plt.ylabel('log volume')
    plt.xlim([-6,4])
    outname = 'p19_pdf_u07_u11_u12.pdf'
    plt.savefig(outname)
    print outname

for frame in  [50,100]: #[0,5,10,20,30,40,50]:
    if 'proj' not in dir() or True:
        #basedir = '/scratch1/dcollins/Paper19/u07-64-driving-nograv'; name='u07'
        #basedir = '/scratch1/dcollins/Paper19/u11_u07_newtracer'; name = 'u11'
        basedir = '/scratch1/dcollins/Paper19/u12_montecarlo'; name = 'u12'
        fname = '%s/DD%04d/data%04d'%(basedir,frame,frame)

        ds = yt.load(fname)


    if 1:
        ad = ds.all_data()
        ad.set_field_parameter('particle_norm',8*ds.index.grids[0]['particle_mass'][0]/(ds.index.grids[0].dds.prod()))
        tn = ad['tracer_number_density']
        phist(tn)
    if 1:

        den_min = 0.1 # min([0.5,ad['density'].min(), ad['tracer_number_density'].min()])
        den_max = 10 #max([5.0,ad['density'].max(), ad['tracer_number_density'].max()])
        phase = yt.create_profile(ad,bin_fields=['density','tracer_number_density'], fields=['cell_mass'],
                                  weight_field=None,
        extrema={'density':[den_min, den_max], 'tracer_number_density':[den_min,den_max]}, n_bins=[128,128])
        pp = yt.PhasePlot.from_profile(phase)
        pp.set_xlabel('Volume Density')
        pp.set_ylabel('CIC tracer number density')
        print pp.save('p19_%s_test_phase_n%04d.pdf'%(name,frame))
#       phase = yt.create_profile(ad,bin_fields=['density'], fields=['cell_mass'],
#                                 weight_field=None,
#       extrema={'density':[den_min, den_max], 'tracer_number_density':[den_min,den_max]}) #, n_bins=[128,128])


    if 1:
        proj =ds.proj('density',0,center=[0.5]*3)
        pw = proj.to_pw(center=[0.5]*3)
        pw.set_cmap('density','gray')
        pw.annotate_particles(1.0,col='r',stride=16)
        print pw.save('p19_%s_n%04d.pdf'%(name,frame))
        #frb = proj.to_frb(1,[64,64])

    if 0:
        plt.clf()
        plt.imshow(np.log10(frb['density']),cmap='gray',origin='lower',interpolation='nearest')
        plt.savefig('p19_frb_test_density_%04d'%(frame))
        plt.clf()
        plt.imshow(np.log10(frb[('deposit','all_cic')]),cmap='gray',origin='lower',interpolation='nearest')
        plt.savefig('p19_%s_frb_test_cic_%04d'%(name,frame))

    #for n,x in enumerate(np.arange(0.,1.,1./32)):
    for n,x in enumerate([]):
        slc = ds.slice(0,x)
        frb=slc.to_frb(1,[64]*2)
        for j,field in enumerate([('deposit','all_cic'), 'density']):
            plt.clf()
            plot=plt.imshow(np.log10(frb[field]),origin='lower',interpolation='nearest')
            plt.colorbar(plot)
            outname = 'p19_test_frb_%s_n%04d_sl%04d'%(['trac','dens'][j],frame,n)
            print outname
            plt.savefig(outname)
