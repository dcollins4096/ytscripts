from go import *
if '/home/dcollins/repos/p49c/p49_eigenmodes' not in sys.path:
    sys.path.append('/home/dcollins/repos/p49c/p49_eigenmodes')
if './p49c' not in sys.path:
    sys.path.append('./p49c')

import p49_eigen
reload(p49_eigen)
import p49_plot_tools
reload(p49_plot_tools)
import p49_QU2EB
import turb_quan
reload(turb_quan)
import plots
reload(plots)

from p49_print_tools import *
def dumb(Q,k):
    x = np.arange(Q.size)
    N = Q.size
    dumbest =  Q*np.exp(-2*np.pi*1j*k*x/N)
    plt.clf()
    plt.plot(1.*x/N,dumbest)
    ndumb = len(glob.glob('p49_dumb/*png'))
    tots=np.sum(dumbest)
    plt.title("%0.2e + i %0.2e"%(tots.real,tots.imag))

    plt.savefig('p49_dumb/dumb_%02d.png'%ndumb)


    return tots
def symmetric(v):
    #for real fft, Ak = A(-k)^*; the negative-phase is the conjugate
    #(that way when you sum, the imaginary parts cancel.)
    #This lets us take an arbitrary K-space signal and ensure it's inverse-fft is real.
    s2=np.zeros_like(v)
    Nx=v.size
    s2[1:Nx/2] = v[1:Nx/2]
    s2[Nx:Nx/2:-1] = v[1:Nx/2].conj()
    s2[0]=v[0]
    return s2


def pp(Q,msg=''):
    print('===')
    print("%s R %s"%(msg,str(Q.real)))
    print("%s I %s"%(msg,str(Q.imag)))

def Nmodes(A=[],B=[],K=[],Theta=[],linthreshy=1e-8, real_fft=True):
    output={}
    size=32
    nmodes=len(A)
    twopi=np.pi*2
    x = np.mgrid[0:1:1./size]
    consts=np.zeros([nmodes])
    linears=np.zeros([nmodes,size])
    sums=np.zeros([nmodes,size])
    for n in range(nmodes):
        consts[n]=A[n]
        linears[n,:]=B[n]*(np.exp(1j*(2*np.pi*K[n]*x+Theta[n])))
        sums[n,:] = consts[n]+linears[n,:]

    Q_flat = np.prod(sums,axis=0)
    plt.clf()
    for n in range(nmodes):
        plt.plot(x,sums[n,:],label=n)
    plt.plot(x,Q_flat,label='prod')
    plt.legend(loc=0)

    plt.savefig('p49c_plots/fftN_1d_test.png')
    #Q_flat = symmetric(Q_flat)



    if real_fft:
        actual_fft=np.fft.rfft(Q_flat)
    else:
        actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(x.size)

    fig,ax = plt.subplots(1,2,figsize=(8,4))
    ax0=ax[0]
    ax1=ax[1]

    ax0.plot(qfft.real,label='r')
    ax0.plot(qfft.imag,label='i')
    ax0.set_yscale('symlog', linthreshy=linthreshy)
    nonzero_k= nz(qfft)
    nonzero_v = nonzero(qfft)

    print('nonzero k %s'%str(nonzero_k))
    pp(nonzero_v,msg='nonzero values')
    import amp_tools
    reload(amp_tools)
    values={}
    values['A'] = nar(A)
    values['B'] = nar(B)
    values['K'] = nar(K)
    values['Theta']=nar(Theta)

    expect = amp_tools.modes_2(**values)
    output.update(values)

    labs = expect.pop('labs')
    output['labs']=labs
    output['full']=expect.pop('full_history')
    yv = np.logspace(-5,1,len(labs))
    for nmode,l in enumerate(labs):
        my_y=yv[nmode]
        ax0.text( 10,my_y,l)
        ax0.plot(  [labs[l],10], [my_y,my_y],c='k')
    expect_k = nar(sorted(list(expect.keys())))
    expect_v = nar([ expect[ key] for key in expect_k])
    ax0.scatter(expect_k,expect_v.real, label='xr',c='c')
    ax0.scatter(expect_k,expect_v.imag, label='xi',c='r')

    output.update({'Q_flat':Q_flat,'actual_fft':actual_fft,'qfft':qfft})
    output['expect']=expect
    output['nonzero_v']=nonzero_v
    output['nonzero_k']=nonzero_k
    output['expect_v']=expect_v
    output['expect_k']=expect_k
    output['B']=B

    def expected_k_matcher(nz_k,nz_v,ex_k,ex_v):
        out_v = np.zeros_like(nz_v)
        for n,k in enumerate(nz_k):
            if k in ex_k:
                out_v[n] = ex_v[ np.where(ex_k == k)]
        return out_v

    if expect_k.size != nonzero_k[0].size:
        print("K sizes don't match")
        ok_expect = expected_k_matcher(nonzero_k[0],nonzero_v,expect_k,expect_v)
        ax1.scatter(nonzero_k, (ok_expect-nonzero_v).imag,label='imag',c='c')
     	ax1.scatter(nonzero_k, (ok_expect-nonzero_v).real,label='real',c='y')
        ax1.set_title("K error. L1 %0.2e"%(np.abs(ok_expect-nonzero_v).sum()))
    else:
        if np.abs( expect_k-nonzero_k).sum() > 1e-7:
            print("k values don't match!!!")
        else:
            print( "Total difference: %0.2e"%( np.abs( expect_v-nonzero_v).sum()))
            ax1.scatter(expect_k, (expect_v-nonzero_v).imag,label='imag',c='r')
            ax1.scatter(expect_k, (expect_v-nonzero_v).real,label='real',c='b')
            pp(expect_v,msg='expected values')
            pp(expect_v-nonzero_v,msg='difference')
            ax1.set_title("KOK. L1 %0.2e"%(np.abs(expect_v-nonzero_v).sum()))
            if np.abs(expect_k-nonzero_k).sum() > 0:
                print("VERY TERRIBLE ERROR K ERROR")
                ax1.set_title("SHIT K ERROR")


    #ax0.legend(loc=0)
    ax1.set_yscale('symlog', linthreshy=1e-8)

    fig.savefig('p49c_plots/fftN_1d_test2.png')
    plt.close(fig)
    return output
vals_small={   "A0" : 1, "A1" : np.sqrt(2.e-6), "KA":1, "Atheta" : 0,
   "B0" : 1, "B1" : np.sqrt(2.e-6), "KB":1, "Btheta" : 0, 
   "C0" : 1, "C1" : 0, "KC":5, "Ctheta":1.2}

vals_0={"A0" : 8, "A1" : 2, "KA":1, "Atheta" : 0,
        "B0" : 2, "B1" : 4, "KB":1, "Btheta" : 0, 
        "C0" : 1, "C1" : 0, "KC":5, "Ctheta":1.2}
vals_n = {'A':nar([1,1,1]),'B':nar([1e-2,1e-2,1e-2]),'K':nar([1,1,5]),'Theta':nar([0,1.2,1.2])}
vals_n = {'A':nar([0,0,1]),'B':nar([1e-2,1e-2,0000]),'K':nar([1,1,0]),'Theta':nar([0,0,0])}
vn = {'A':nar([1,1.4]),'B':nar([1,1e-2]),'K':nar([1,1]),'Theta':nar([0,np.pi/1.2])}
vn = {'A':nar([1,1.4]),'B':nar([1,1e-2]),'K':nar([1,3]),'Theta':nar([0,np.pi/1.2])}
vn = {'A':nar([1,0,1]),'B':nar([1,1e-2,1e-2]),'K':nar([1,1,3]),'Theta':nar([0.1,0.3,0.5])}
#vn = {'A':nar([1,0,1,0]),'B':nar([1,1e-2,1e-2,1e-3]),'K':nar([1,1,3,4]),'Theta':nar([0.1,0.3,0.5,0])}
#vn = {'A':nar([1,1,1,0]),'B':nar([1,0   ,0    ,1e-3]),'K':nar([1,1,3,4]),'Theta':nar([0.1,0.3,0.5,0])}
vn = {'A':nar([0,0,0,1]),'B':nar([1,1   ,1    ,1]),'K':nar([1,2,3,4]),'Theta':nar([0.2,0.1,np.pi/2,0])}
vn = {'A':nar([0,0,0,0,0]),'B':nar([1,1,1,1,1]),'K':nar([1,1,1,1,1]),'Theta':nar([0,0,0,0,0])}#nar([0,0.2,0.1,np.pi/2,0])}
vn = {'A':nar([1,1,1,1,1]),'B':nar([1,1,1,1,1]),'K':nar([1,1,1,1,1]),'Theta':nar([0.2,0.1,np.pi/2,0,1.3,  4])}#nar([0,0.2,0.1,np.pi/2,0])}

vals_1={"A0" : vn['A'][0], "A1" :vn['B'][0] , "KA":vn['K'][0], "Atheta" :vn['Theta'][0],
        "B0" : vn['A'][1], "B1" :vn['B'][1] , "KB":vn['K'][1], "Btheta" :vn['Theta'][1],
        "C0" : 1         , "C1" :0          , "KC":0         , "Ctheta" :0             }
if len(vn['A'] ) > 2:
    vals_1.update({"C0" : vn['A'][2], "C1" :vn['B'][2] , "KC":vn['K'][2], "Ctheta" :vn['Theta'][2]})
oneN=Nmodes(linthreshy=0.1,**vn)
#oned=stuffs(**vals_1)
