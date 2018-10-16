
import p44_random_tools as RT
reload(RT)
if 0:
    plt.clf()
    t = np.arange(0,np.pi*2,0.01)
    s = np.sin(t)
#plt.plot(t,s)
#plt.savefig('p44_t_h.png')
    plt.hist(s, histtype='step',bins=100)
    outname = 'p44_t_h_hist.png'
    plt.savefig(outname)
    print(outname)

Nx = 100
k  = np.arange(Nx,dtype='complex')
x = np.arange(0,1,1./Nx)
if 1:

    Ak = np.zeros(Nx,dtype='complex')
    odds = k.real%2 == 1
    # -Nx keeps the normalization right, so F(x) = 1,0
    # 1/np.pi is for the actual series.
    # 1j makes it a sign series.
    # Ak[0]=50 keeps the zero-point right (otherwise it's +- 1/2)
   
    Ak[odds] = -Nx/np.pi/k[odds]*1j #1/2 + 2./np.pi/k[odds]
    Ak[0] += 50
    Ak = RT.symmetric(Ak)
    ax = np.fft.ifft(Ak)
    tp = np.abs(ax)
    plt.clf()
    dumb_plt(plt,None,tp,'x','square','p44_square_test.pdf')

    #verify with an actual square wave
    rect_fft = np.fft.fft(x<1/2.)
    plt.clf()
    dumb_plt(plt,None,np.abs(rect_fft), 'k','|Ak|','p44_square_test_k.png')
    dumb_plt(plt,None,np.abs(Ak), 'k','|Ak|','p44_square_test_k.png',c='r')

    #Just two modes, to ensure the normalization is right.
    two = F_n(2)
    plt.clf()
    dumb_plt(plt,x,two,'x','two','p44_square_test_F_2.png')
    two_hat = np.fft.fft(two)
    plt.clf()
    dumb_plt(plt,None,two_hat.real,'k','Two hat','p44_square_test_F_2_hat.png',c='r')
    dumb_plt(plt,None,two_hat.imag,'k','Two hat','p44_square_test_F_2_hat.png',c='g')
    #plt.plot(plt,None,two_hat.imag,'k','Two hat','p44_square_test_F_2_hat.png',c='g')

if 0:
    sin_n = lambda n: np.sin(2 * np.pi * n * x)
    F_n = lambda n: (1/2. + np.sum([2./np.pi/(2*i+1) * np.sin(2 * np.pi * (2*i+1) * x) for i in range(n+1)], axis=0))
    plt.clf()
    for i in [0, 1, 2, 3, 5, 10]:
            plt.plot(x, F_n(i))
            plt.ylim(-1.1, 2.1)
    dumb_plt(plt,x, x<1/2.,'x','square','p44_square_test_2.png')
    plt.clf()
    plt.clf()
    sum_i_j = lambda i,j: np.sum([sin_n(k) for k in range(i, j+1)], axis=0)
    for i in [1, 2, 3, 4, 5, 10, 20, 40]:
            plt.plot(x, sum_i_j(1, i))
    plt.savefig('p44_square_test_fail.png')

    

    
