import amp_tools
reload(amp_tools)
if 0:
    size=32
    twopi=np.pi*2
    x,y,z = np.mgrid[0:1:1./size, 0:1:1./size,0:1:1./size]
    k_unit_int = nar([1,7,1])
    k_unit_Q = np.array(k_unit_int)*twopi
    ampl=0.1
    Q_flat =  ampl* (np.exp(1j*(k_unit_Q[0]*x+k_unit_Q[1]*y+k_unit_Q[2]*z))).real
    #k_unit_Q = np.array([1,1])*twopi
    #Q_flat +=   0.1*(np.exp(1j*(k_unit_Q[0]*x+k_unit_Q[1]*y))).imag
    qfft = np.fft.rfftn(Q_flat)
    print("Nonzero at %s"%str(k_unit_int))
    print("should match %s"%str(nz(qfft)))
    print("amplitude %s"%str(ampl))
    nonzeros = nonzero(qfft)/x.size
    print("should match %s"%str(nonzeros.real))
    print("should be zero %s"%str(nonzeros.imag))
    print("somehow off by 0.5")

    qfft2 =  np.fft.rfftn(Q_flat**2)
    print("should also be twice the unit vectors, and zero %s"%str(2*k_unit_int))
    print("match %s"%str(nz(qfft2)))

if 0:
    import scipy.signal
    shape=5
    #a1 = np.zeros(shape)
    #a1[2] = 1
    ##a1 = np.arange(shape)
    #a2 = np.arange(shape)
    a1=np.arange(shape)
    a2=np.arange(shape)*10
    a1=np.zeros(shape)
    a2=np.zeros(shape)
    a1[1]=10
    a2[1]=20

    b1 = scipy.signal.convolve(a1,a2,mode='full')
    b2 = scipy.signal.convolve(a1,a2,mode='valid')
    b3 = scipy.signal.convolve(a1,a2,mode='same')
    b4 = scipy.signal.convolve(a1,a2,mode='same',method='fft')
    print("a1 %s"%str(a1))
    print("a2 %s"%str(a2))
    #print("b1 %s"%str(b1))
    #print("b2 %s"%str(b2))
    #print("b3 %s"%str(b3))
    #print("b4 "+ " %0.2f"*shape%tuple(b4))
    pnz(b4)

    other_thing = np.zeros(shape)
    if 0:
        for n in range(shape):
            print("===",n)
            a2hat =  np.roll(a2,n)[::-1]
            prod = a1*(a2hat)
            print("a1  = %s"%str(a1))
            print("a2- = %s"%str(a2hat))
            print("pro = %s"%str(prod))
            print("tot = %s"%str(other_thing))
            other_thing[n] = (prod).sum()
        print("===")
        print(other_thing)



if 0:
    #this works for the projection.
    #Doing this for the vector fields is very difficult.
    print('=== density')
    ok=nz(stuff['s'].kstuff['ampl']) 
    amp=nonzero(stuff['s'].kstuff['ampl']) 
    d_pert=stuff['s'].rot['d'][ok]*amp
    dhat = nonzero(stuff['hats']['dhat'])/32.**3
    d0=dhat[0]
    d1=dhat[1]
    print("Density: this should be one",np.abs(dhat[1]/(d_pert)))
    print("values",np.abs(dhat/32**3))

if 0:
    #this works for the projection.
    print('=== h hori')
    ok=nz(stuff['s'].kstuff['ampl']) 
    amp=nonzero(stuff['s'].kstuff['ampl']) 

    h_pert=stuff['s'].rot['hy'][ok]
    hhat = nonzero(stuff['hats']['hhhat'])
    h_expected = h_pert*amp*32**3
    #print("H works",nonzero(np.fft.rfftn(stuff['Hh']))/h_pert/32**3)
    print( "Should be 1, hh: %0.2f"% (hhat[1]/(h_expected)).real)
    print(np.abs(hhat/32**3))

    print('=== h vert')
    v_pert=stuff['s'].rot['hz'][ok]
    vhat = nonzero(stuff['hats']['hvhat'])
    v_expected = v_pert*amp*32**3
    #print("H works",nonzero(np.fft.rfftn(stuff['Hh']))/h_value/32**3)
    print( "Should be 1 hv: %0.2f"% (vhat[1]/(v_expected)).real)
    pnz(np.abs(vhat/32**3))

    print('=== b^2')
    bx_pert=stuff['s'].rot['hx'][ok]*amp
    by_pert=stuff['s'].rot['hy'][ok]*amp
    bz_pert=stuff['s'].rot['hz'][ok]*amp

    bx_hat = np.fft.rfftn(stuff['s'].cubes['hx'])
    bx_nonzero = nonzero(bx_hat)/32.**3
    bx0=bx_nonzero[0]
    bx0_test=bx0/stuff['s'].b0[0]
    bx1=bx_nonzero[1]
    bx1_test = bx1/bx_pert

    by_hat = np.fft.rfftn(stuff['s'].cubes['hy'])
    by_nonzero = nonzero(by_hat)/32.**3
    by0=by_nonzero[0]
    by0_test=by0/stuff['s'].b0[1]
    by1=by_nonzero[1]
    by1_test = by1/by_pert

    bz_hat = np.fft.rfftn(stuff['s'].cubes['hz'])
    bz_nonzero = nonzero(bz_hat)/32.**3
    bz0=bz_nonzero[0]
    bz0_test=bz0/stuff['s'].b0[2]
    bz1=bz_nonzero[1]
    bz1_test = bz1/bz_pert
    b0_test = nar([bx0_test,by0_test,bz0_test])
    b1_test = nar([bx1_test,by1_test,bz1_test])

#   #print("v1 %0.2e v2 %0.2e"%(b2_value, stuff['H2'].max()))
#   b2_expected = b2_value*amp*32**3
#   b2_nonzero = nonzero(stuff['hats']['b2ha'])/32.**3

#   b20_test = b2_nonzero[0]/(bx0*bx0+by0*by0+bz0*bz0)
#   b21_test = b2_nonzero[1]/(2*(bx0*bx1+by0*by1+bz0*bz1))
#   b22_test = b2_nonzero[2]/((bx1*bx1+by1*by1+bz1*bz1))
#   b_total_test = nar([b20_test,b21_test,b22_test])

if 0:
    """this is correct."""
    bx2_xspace=stuff['s'].cubes['hx']*stuff['s'].cubes['hx']
    bx2_hat = np.fft.rfftn(bx2_xspace)
                          #stuff['s'].cubes['hy']*stuff['s'].cubes['hy']+
                          #stuff['s'].cubes['hz']*stuff['s'].cubes['hz'])
    bx2_nonzero = nonzero(bx2_hat)/bx2_xspace.size
    bx2_0_expect = (bx0*bx0+2*bx1*bx1)
    bx2_1_expect = (2*(bx0*bx1))
    bx2_2_expect = ((bx1*bx1))
    bx20_test = bx2_nonzero[0]/bx2_0_expect
    bx21_test = bx2_nonzero[1]/bx2_1_expect
    bx22_test = bx2_nonzero[2]/bx2_2_expect
    bx2_all_expect = nar([bx2_0_expect,bx2_1_expect,bx2_2_expect])
    bx2_all_test = nar([bx20_test,bx21_test,bx22_test])
    print("FT( B^2) %s"%str(bx2_nonzero.real))
    print("FT( B^2) e %s"%str(np.abs(bx2_all_expect)))
    print("FT( B^2) t %s"%str(np.abs(1-bx2_all_test)))

if 0:
    print("QQQ")
    qf = nonzero(stuff['Qh'])/stuff['Delta'][0]**2/32.**2 
    qhatb = np.fft.rfftn(stuff['Q'])/32.**2 
    diff = 5 * np.pi / 180/32.
    qhatbA = nonzero(qhatb)
    #qf = nonzero(stuff['Q']) #/32.**3 don't need this, already done.
    print('Qf',qf)
    q0 = d0*(by0*by0-bz0*bz0)
    print('Qf0 test', q0, qf[0])
    turb_quan.plotter2([ np.fft.irfftn( stuff['Qh']), stuff['Q']],'p49c_plots/test_Qh_Q.png',
                       norm='ind')


#some dumb tests0
def moof(arr):
    return( np.fft.rfftn(arr))
def moo(arr):
    return( nonzero(np.fft.rfftn(arr))/arr.size)
def monz(arr):
    return( nz(np.fft.rfftn(arr)))
def pp(Q,msg=''):
    print('===')
    print("%s R %s"%(msg,str(Q.real)))
    print("%s I %s"%(msg,str(Q.imag)))
    
c=stuff['s'].cubes
if 'hatd' not in dir():
    """this works"""
    c=stuff['s'].cubes
    ql2=c['d']*c['d'] #*( c['hy']**2-c['hz']**2)
    hatl2 =  moo(ql2)
    hatd = moo( c['d'])
    hathy = moo( c['hy'])
    hathx = moo( c['hx'])
    pp(1- hatl2[0]/( hatd[0]**2 + 2*hatd[1]**2), 'dd relerr k=0')
    pp(1-hatl2[1]/( hatd[1]*hatd[0]*2), 'dd relerr k=1')
    pp(1-hatl2[2]/( hatd[1]**2), 'dd relerr k=2')

if 0:
    """this works"""
    a2=c['d']#*c['d'] #*c['hy'] #*( c['hy']**2-c['hz']**2)
    moo_a2_complex = moo(a2)
    moo_a2=moo_a2_complex.real
    KA = 1
    KB = 1
    KC = 1
    A0 = hatd[0]; A1 = 2*hatd[1]
    B0 = 1; B1=0; #hatd[0]; B1 = hatd[1]
    #B0 = hathy[0]; B1 = 2*hathy[1]
    C0 = 1; C1 = 0
    expect_d1 = amp_tools.modes_1(KA,KB,KC,A0,A1,B0,B1,C0,C1)
    labs = expect_d1.pop('labs')
    expect_k = nar(sorted(list(expect_d1.keys())))
    expect_v = nar([ expect_d1[ key] for key in expect_k])
    for k in range(len(moo_a2)):
        print("==",k)
        print('x',expect_v[k])
        print('a',moo_a2[k])
        print('err', 1-expect_v[k]/moo_a2[k])

if 0:
    """works"""
    d2=c['d']*c['d'] #*c['hy'] #*( c['hy']**2-c['hz']**2)
    moo_d2_complex = moo(d2)
    moo_d2=moo_d2_complex.real
    values={}
    values['KA']=1
    values['KB']=1
    values['KC']=1
    values['A0']=hatd[0]
    values['A1']=2.*hatd[1]
    values['B0']=hatd[0]
    values['B1']=2.*hatd[1]
    values['C0']=1
    values['C1']=0
    values['Atheta']=0
    values['Btheta']=0
    values['Ctheta']=0
    expect_d2 = amp_tools.modes_1(**values)
    labs = expect_d2.pop('labs')
    full = expect_d2.pop('full_history')
    ex_d2_k = nar(sorted(list(expect_d2.keys())))
    ex_d2_v = nar([ expect_d2[ key] for key in ex_d2_k])
    def dfer(xpct):
        for k in range(len(moo_d2)):
            print("==",k)
            print('x',xpct[k])
            print('a',moo_d2[k])
            print('err', 1-xpct[k]/moo_d2[k])
    dfer(ex_d2_v)

if 0:
    """this works"""
    a2=c['d']*c['hy'] #*( c['hy']**2-c['hz']**2)
    moo_a2_complex = moo(a2)
    moo_a2=moo_a2_complex.real
    pp(1-moo_a2[0]/( hatd[0]*hathy[0] + 2*hatd[1]*hathy[1]  ), "dhy relerr k=0")
    pp(1-moo_a2[1]/( hatd[0]*hathy[1] + hatd[1]*hathy[0] ), "dhy relerr k=1")
    pp(1-moo_a2[2]/( hatd[1]*hathy[1]), "dhy relerr k=2")

if 0:
    """works."""
    dhy=c['hy']*c['d'] #*c['hy'] #*( c['hy']**2-c['hz']**2)
    moo_dhy_complex = moo(dhy)
    moo_dhy=moo_dhy_complex.real
    values={}
    values['KA']=1
    values['KB']=1
    values['KC']=1
    values['A0']=hatd[0]
    values['A1']=2.*hatd[1]
    values['B0']=hathy[0]
    values['B1']=2.*hathy[1]
    values['C0']=1
    values['C1']=0
    values['Atheta']=0
    values['Btheta']=0
    values['Ctheta']=0
    expect_dhy = amp_tools.modes_1(**values)
    labs = expect_dhy.pop('labs')
    full = expect_dhy.pop('full_history')
    ex_dhy_k = nar(sorted(list(expect_dhy.keys())))
    ex_dhy_v = nar([ expect_dhy[ key] for key in ex_dhy_k])
    def dfer(xpct,actual):
        for k in range(len(actual)):
            print("==",k)
            print('x',xpct[k])
            print('a',actual[k])
            print('err', 1-xpct[k]/actual[k])
    dfer(ex_dhy_v, moo_dhy)
    print("total error", np.abs(ex_dhy_v-moo_dhy).sum())

if 1:
    """in progress."""
    dhyhy=c['hy']*c['d']*c['hy'] #*( c['hy']**2-c['hz']**2)
    moo_dhyhy_complex = moo(dhyhy)
    moo_dhyhy=moo_dhyhy_complex.real
    values={}
    values['KA']=1
    values['KB']=1
    values['KC']=1
    values['A0']=hatd[0]
    values['A1']=2.*hatd[1]
    values['B0']=hathy[0]
    values['B1']=2.*hathy[1]
    values['C0']=hathy[0]
    values['C1']=2.*hathy[1]
    values['Atheta']=0
    values['Btheta']=0
    values['Ctheta']=0
    expect_dhyhy = amp_tools.modes_1(**values)
    labs = expect_dhyhy.pop('labs')
    full = expect_dhyhy.pop('full_history')
    ex_dhyhy_k = nar(sorted(list(expect_dhyhy.keys())))
    ex_dhyhy_v = nar([ expect_dhyhy[ key] for key in ex_dhyhy_k])
    def dfer(xpct,actual):
        for k in range(len(actual)):
            print("==",k)
            print('x',xpct[k])
            print('a',actual[k])
            print('err', 1-xpct[k]/actual[k])
    dfer(ex_dhyhy_v, moo_dhyhy)
    print("total error", np.abs(ex_dhy_v-moo_dhy).sum())

if 0:
    """this doesn't"""
    a2=c['d']*c['hy'] #*c['hy'] #*( c['hy']**2-c['hz']**2)
    moo_a2_complex = moo(a2)
    moo_a2=moo_a2_complex.real
    KA = 1
    KB = 1
    KC = 1
    A0 = hatd[0]; A1 = 2*hatd[1]
    B0 = hathy[0]; B1 = 2*hathy[1]
    C0 = 1; C1 = 0
    expect = amp_tools.modes_1(KA,KB,KC,A0,A1,B0,B1,C0,C1)
    labs = expect.pop('labs')
    expect_k = nar(sorted(list(expect.keys())))
    expect_v = nar([ expect[ key] for key in expect_k])
    
if 0:
    a2=c['d']*c['hy'] #*( c['hy']**2-c['hz']**2)
    moo_a2_complex = moo(a2)
    moo_a2=moo_a2_complex.real
    pp(1-moo_a2[0]/( hatd[0]*hathy[0] + 2*hatd[1]*hathy[1]  ), "dhy relerr k=0")
    pp(1-moo_a2[1]/( hatd[0]*hathy[1] + hatd[1]*hathy[0] ), "dhy relerr k=1")
    pp(1-moo_a2[2]/( hatd[1]*hathy[1]), "dhy relerr k=2")


    pp(1-moo_a2[0]/(expect_v[0]), "butts k=0")
    pp(1-moo_a2[1]/(expect_v[1]), "butts k=1")
    pp(1-moo_a2[2]/(expect_v[2]), "butts k=2")


