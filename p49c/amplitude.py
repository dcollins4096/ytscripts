
if 1:
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

    thing = scipy.signal.convolve(qfft,qfft,method='fft')
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



"""
if 1:
    #this works for the projection.
    #Doing this for the vector fields is very difficult.
    ok=nz(stuff['s'].kstuff['ampl']) 
    amp=nonzero(stuff['s'].kstuff['ampl']) 
    d_value=stuff['s'].rot['d'][ok]
    dhat = nonzero(stuff['hats']['dhat'])
    print("Density: this should be one",dhat[1]/(d_value*amp*32.**3))
if 1:
    #this works for the projection.
    ok=nz(stuff['s'].kstuff['ampl']) 
    amp=nonzero(stuff['s'].kstuff['ampl']) 

    h_value=stuff['s'].rot['hy'][ok]
    hhat = nonzero(stuff['hats']['hhhat'])
    h_expected = h_value*amp*32**3
    #print("H works",nonzero(np.fft.rfftn(stuff['Hh']))/h_value/32**3)
    print( "Should be 1, hh: %0.2f"% (hhat[1]/(h_expected)).real)

    v_value=stuff['s'].rot['hz'][ok]
    vhat = nonzero(stuff['hats']['hvhat'])
    v_expected = v_value*amp*32**3
    #print("H works",nonzero(np.fft.rfftn(stuff['Hh']))/h_value/32**3)
    print( "Should be 1 hv: %0.2f"% (vhat[1]/(v_expected)).real)

    bx_value=stuff['s'].rot['hx'][ok]
    by_value=stuff['s'].rot['hy'][ok]
    bz_value=stuff['s'].rot['hz'][ok]
    b2_value = (bx_value**2+by_value**2+bz_value**2)
    print("v1 %0.2e v2 %0.2e"%(b2_value, stuff['H2'].max()))
    b2_expected = b2_value*amp*32**3
    b2hat = nonzero(stuff['hats']['b2ha'])
    #print("H works",nonzero(np.fft.rfftn(stuff['Hh']))/h_value/32**3)
    print( "Should be 1 b2: %0.2f"% (b2hat[1]/(b2_expected)).real)

if 1:
    thresh = 1e-11
    flat_mask = np.abs(stuff['hats']['hvhat']) > thresh
    q_flat = stuff['Qh'][flat_mask]
    q_prolly = (dhat * (h_value**2 - v_value**2 )) /(amp*32.**3)
    print(q_flat/q_prolly)
    """
