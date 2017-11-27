if 'ef' not in dir(): execfile('go')

import p49_fields
reload(taxi)
if 'car' not in dir():
    car=taxi.taxi('g18e_tight')
    ds=car.load()

if 0:
    g = ds.index.grids[-1]
    #N2name = 'dave1_Qz_n0-0001_p-1' # 'dave1_Qz_n0-0001_p-1' #"dave4_N2z_n0-0001_p-1"
    N2name = 'Uz_n0-0001_p-1' # 'dave1_Qz_n0-0001_p-1' #"dave4_N2z_n0-0001_p-1"
    n=g['density'].in_units('code_density')
    n2 = g[N2name].in_units('code_density')
    stat(n,'density')
    stat(n2,N2name)

if 0:
    car.fields=['density', N2name] #,'magnetic_field_strength']
    car.callbacks=[]
    car.axis=[2]
    car.plot()

if 'the_plot' not in dir(): # and False:
    car=taxi.taxi('g18e_tight')
    car.operation='CenterSlice'
    #car=taxi.taxi('gdl')
    car.fields=['density'] #, N2name] #,'magnetic_field_strength']
    Qname = "Qz_n0-0001_p-1"
    #car.fields = [Qname]
    #car.operation='CenterSlice'
    car.axis=[2]
    #car.callbacks = ['stokes_angles']
    car.callbacks = []
    car.outname = 'g30_%s_stoke2'%(car.name)
    kludge_my_package={}
    car.plot()
    this_frb = car.the_plot.frb
    the_plot=car.the_plot


if 1:
    ef('QU_callback.py')
    the_plot._callbacks=[]
    #the_plot._callbacks.append(StokesCallback(weight='N', kludge_my_package=kludge_my_package, normalize=True, scale=50))
    the_plot.annotate_stokes_angles2(weight='I', kludge_my_package=kludge_my_package, normalize=False, scale=50)
    the_plot.run_callbacks()
    print the_plot.save('g30_%s_annotate_6.png'%car.name)
    print "did the thing"
if 0:
    plt.clf()
    Uname = "Uz_n0-0001_p-1"
    Qname = "Qz_n0-0001_p-1"
    #Uname = "Uz_n0-0001_p-1"
    #Qname = "Qz_n0-0001_p-1"
    
    N2name ="N2z_n0-0001_p-1"
    #Qname = "Q%s_n0-%04d_p-%d"%(axis_dict[ax], self.n0, self.p)
    #N2name = "N2%s_n0-%04d_p-%d"%(axis_dict[ax], self.n0, self.p)
#length = np.ones_like(theta_stokes)*self.length
    Q = this_frb[Qname]
    U = this_frb[Uname]
    N = this_frb['density']#.in_units('code_density').v
    N2 =this_frb[N2name]
    xx0=0; xx1=800; nx=800
    yy0=0; yy1=800; ny=800

    theta_stokes = 0.5*np.arctan2(U,Q)
    p0 = 0.1
    frac = p0*np.sqrt(Q**2+U**2)/(N - p0*N2)
    frac = N/N.max()
    fv_x = frac*np.cos(theta_stokes)
    fv_y = frac*np.sin(theta_stokes)

    #plt.clf()
    #a=plt.imshow(frac, interpolation='nearest',origin='lower')#, cmap='gray')
    #plt.colorbar(a)
    #plt.savefig('g30b_test_frac.png')

    plt.clf()
    a=plt.imshow(p0*np.sqrt(Q**2+U**2), interpolation='nearest',origin='lower')#, cmap='gray')
    plt.colorbar(a)
    plt.savefig('g30b_test_I.png')

    plt.clf()
    a=plt.imshow(p0*N2, interpolation='nearest',origin='lower')#, cmap='gray')
    plt.colorbar(a)
    plt.savefig('g30b_test_N2.png')

    plt.clf()
    a=plt.imshow(N-p0*N2, interpolation='nearest',origin='lower')#, cmap='gray')
    plt.colorbar(a)
    plt.savefig('g30b_test_N-N2.png')

    plt.clf()
    a=plt.imshow(np.sqrt(Q**2+U**2)/N, interpolation='nearest',origin='lower')#, cmap='gray')
    plt.colorbar(a)
    plt.savefig('g30b_test_IoverN.png')
    
    plt.imshow(np.log10(N), interpolation='nearest',origin='lower')#, cmap='gray')
    X,Y = np.meshgrid(np.linspace(xx0,xx1,nx,endpoint=True),
                      np.linspace(yy0,yy1,ny,endpoint=True))
    sl1 = slice(None,None,10)
    sl = [sl1,sl1]
    stat(X[sl])
    #mp=car.kludge_my_package

    whatzit = np.sqrt( (fv_x[sl]**2 + fv_y[sl]**2).max())
    plt.quiver(X[sl],Y[sl], fv_x[sl], fv_y[sl], headlength=0, headwidth=1,color='r',scale=50./whatzit)
    outname = 'g30b_test_N_stokes_frac4.png'
    plt.savefig(outname); print outname
    plt.clf()
