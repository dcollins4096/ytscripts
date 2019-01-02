
if 0:
    size=32
    twopi=np.pi*2
    x,y = np.mgrid[0:1:1./size, 0:1:1./size]
    k_unit_Q = np.array([1,7])*twopi
    Q_flat =  2+0.1* (np.exp(1j*(k_unit_Q[0]*x+k_unit_Q[1]*y))).imag
    #k_unit_Q = np.array([1,1])*twopi
    #Q_flat +=   0.1*(np.exp(1j*(k_unit_Q[0]*x+k_unit_Q[1]*y))).imag
    qfft = np.abs( np.fft.rfftn(Q_flat))
    qfft2 = np.abs( np.fft.rfftn(Q_flat**2))
    print(nz(qfft))

    print(nonzero(qfft)/1024)
    print(nz(qfft2))
    print(nonzero(qfft2)/1024)
    #pnz(qfft)
    #pnz(qfft2)

    #thing = scipy.signal.convolve(qfft,qfft,method='fft')
    #    'thing' just doubles if qfft is a delta.
    #k_dot_x=k_unit_Q[0]*x+k_unit_Q[1]*y
    #qfft_norm = qfft/np.abs(qfft).max()
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
    d_value=stuff['s'].rot['d'][ok]*amp
    dhat = nonzero(stuff['hats']['dhat'])/32.**3
    d0=dhat[0]
    d1=dhat[1]
    print("Density: this should be one",dhat[1]/(d_value))
    print(np.abs(dhat/32**3))
if 0:
    #this works for the projection.
    print('=== h hori')
    ok=nz(stuff['s'].kstuff['ampl']) 
    amp=nonzero(stuff['s'].kstuff['ampl']) 

    h_value=stuff['s'].rot['hy'][ok]
    hhat = nonzero(stuff['hats']['hhhat'])
    h_expected = h_value*amp*32**3
    #print("H works",nonzero(np.fft.rfftn(stuff['Hh']))/h_value/32**3)
    print( "Should be 1, hh: %0.2f"% (hhat[1]/(h_expected)).real)
    print(np.abs(hhat/32**3))

    print('=== h vert')
    v_value=stuff['s'].rot['hz'][ok]
    vhat = nonzero(stuff['hats']['hvhat'])
    v_expected = v_value*amp*32**3
    #print("H works",nonzero(np.fft.rfftn(stuff['Hh']))/h_value/32**3)
    print( "Should be 1 hv: %0.2f"% (vhat[1]/(v_expected)).real)
    pnz(np.abs(vhat/32**3))

    print('=== b^2')
    bx_value=stuff['s'].rot['hx'][ok]*amp
    by_value=stuff['s'].rot['hy'][ok]*amp
    bz_value=stuff['s'].rot['hz'][ok]*amp

    bx_hat = np.fft.rfftn(stuff['s'].cubes['hx'])
    bx_nonzero = nonzero(bx_hat)/32.**3
    bx0=bx_nonzero[0]
    bx0_test=bx0/stuff['s'].b0[0]
    bx1=bx_nonzero[1]
    bx1_test = bx1/bx_value

    by_hat = np.fft.rfftn(stuff['s'].cubes['hy'])
    by_nonzero = nonzero(by_hat)/32.**3
    by0=by_nonzero[0]
    by0_test=by0/stuff['s'].b0[1]
    by1=by_nonzero[1]
    by1_test = by1/by_value

    bz_hat = np.fft.rfftn(stuff['s'].cubes['hz'])
    bz_nonzero = nonzero(bz_hat)/32.**3
    bz0=bz_nonzero[0]
    bz0_test=bz0/stuff['s'].b0[2]
    bz1=bz_nonzero[1]
    bz1_test = bz1/bz_value

    #b2_value = (bx_value**2+by_value**2+bz_value**2)
    #print("v1 %0.2e v2 %0.2e"%(b2_value, stuff['H2'].max()))
    #b2_expected = b2_value*amp*32**3
    b2hat = nonzero(stuff['hats']['b2ha'])/32.**3
#   b20_test = b2hat[0]/(bx0*bx0+by0*by0+bz0*bz0)

#   b21_test = b2hat[1]/(2*(bx0*bx1+by0*by1+bz0*bz1))
#   b22_test = b2hat[2]/((bx1*bx1+by1*by1+bz1*bz1))

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
if 0:
    def moo(arr):
        return( nonzero(np.fft.rfftn(arr))/32.**3)
    c=stuff['s'].cubes
    ql2=c['d']*c['d'] #*( c['hy']**2-c['hz']**2)
    hatl2 =  moo(ql2)
    hatd = moo( c['d'])
    hathy = moo( c['hy'])
    print( hatl2[1]/( hatd[1]*hatd[0]*2))
    print( hatl2[2]/( hatd[1]**2))
    a2=c['d']*c['hy'] #*( c['hy']**2-c['hz']**2)
    moo_a2 = moo(a2).real
    print(1-moo_a2[0]/( hatd[0]*hathy[0] ))
    print(moo_a2[1]/( hatd[0]*hathy[1] + hatd[1]*hathy[0] ))
    print(moo_a2[2]/( hatd[1]*hathy[1]))

if 0:
    def moo(arr):
        return( nonzero(np.fft.rfftn(arr))/32.**3)
    c=stuff['s'].cubes
    hd = moo( c['d'])
    hy = moo( c['hy'])
    a2=c['d']*c['hy']*c['hy'] #*( c['hy']**2-c['hz']**2)
    moo_a2 = moo(a2).real
    print(moo_a2[0]/( hd[0]*hy[0]*hy[0] ))
    print(moo_a2[1]/( hd[0]*hy[0]*hy[1] + 
                      hd[1]*hy[0]*hy[0]+
                      hd[0]*hy[0]*hy[1]))
    print(moo_a2[2]/( hd[0]*hy[1]*hy[1] + 
                      hd[1]*hy[0]*hy[1]+
                      hd[1]*hy[1]*hy[0]))
    print(moo_a2[3]/( hd[1]*hy[1]*hy[1] ))
    
if 0:
    def moo(arr):
        return( nonzero(np.fft.rfftn(arr))/32.**3)
    c=stuff['s'].cubes
    hatd = moo( c['d'])
    hz = moo( c['hz'])
    a3=c['d']*c['hz']*c['hz'] #*( c['hy']**2-c['hz']**2)
    moo_a3 = moo(a3).real
    print(moo_a3[0]/( hatd[0]*hz[0]*hz[0] ))
    print(moo_a3[1]/( hatd[0]*hz[0]*hz[1] + 
                      hatd[1]*hz[0]*hz[0]+
                      hatd[0]*hz[0]*hz[1]))
    print(moo_a3[2]/( hatd[0]*hz[1]*hz[1] + 
                      hatd[1]*hz[0]*hz[1]+
                      hatd[1]*hz[1]*hz[0]))
    print(moo_a3[3]/( hatd[1]*hz[1]*hz[1] ))

if 0:
    a4 = a2-a3
    moo_a4 = moo(a4)
    a5 = c['d']*(c['hy']**2-c['hz']**2)
    moo_a5=moo(a5)
    print(moo_a2)
    print(moo_a3)
    print(moo_a4.real-(moo_a2-moo_a3))
                     
if 0:
    def moo(arr):
        return( nonzero(np.fft.rfftn(arr))/32.**3)
    c=stuff['s'].cubes
    hatd = moo( c['d'])
    hathy = moo( c['hy'])
    hathz = moo( c['hz'])
    a2=c['d']*(c['hy']*c['hy'] -c['hz']*c['hz'])
    moo_a2 = moo(a2).real
    guess = np.zeros(moo_a2.size)
    guess[0]= hatd[0]*hathy[0]*hathy[0] - hatd[0]*hathz[0]*hathz[0]
    print(moo_a2[0]/guess[0])
#   print(moo_a2[1]/( hatd[0]*hathy[0]*hathy[1] + 
#                     hatd[1]*hathy[0]*hathy[0]+
#                     hatd[0]*hathy[0]*hathy[1]))
#   print(moo_a2[2]/( hatd[0]*hathy[1]*hathy[1] + 
#                     hatd[1]*hathy[0]*hathy[1]+
#                     hatd[1]*hathy[1]*hathy[0]))
#   print(moo_a2[3]/( hatd[1]*hathy[1]*hathy[1] ))
#                    
