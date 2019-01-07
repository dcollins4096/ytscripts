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

def pp(Q,msg=''):
    print('===')
    print("%s R %s"%(msg,str(Q.real)))
    print("%s I %s"%(msg,str(Q.imag)))

def stuffs(KA=0,KB=0,KC=0,A0=0,A1=0,B0=0,B1=0,C0=0,C1=0, Atheta=0,Btheta=0,Ctheta=0):
    output={}
    size=32
    twopi=np.pi*2
    x = np.mgrid[0:1:1./size]
    Aa = A0
    Ab = A1* (np.exp(1j*(2*np.pi*KA*x+Atheta))).real
    A =  Aa+Ab

    Ba = B0
    Bb = B1* (np.exp(1j*(2*np.pi*KB*x+Btheta))).real
    B =  Ba+Bb

    Ca = C0
    Cb = C1* (np.exp(1j*(2*np.pi*KC*x+Ctheta))).real
    C =  Ca+Cb
    plt.clf()
    plt.plot(x,A,label='A')
    plt.plot(x,B,label='B')
    #plt.plot(x,C,label='C')
    plt.plot(x,A*B*C,label='A*B*C')
    plt.legend(loc=0)

    plt.savefig('p49c_plots/fft_1d_test.png')
    Q_flat = A*B*C
    actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(x.size)

    fig,ax = plt.subplots(1,2,figsize=(8,4))
    ax0=ax[0]
    ax1=ax[1]

    ax0.plot(qfft.real,label='r')
    ax0.plot(qfft.imag,label='i')
    ax0.set_yscale('symlog', linthreshy=1e-8)
    nonzero_k= nz(qfft)
    nonzero_v = nonzero(qfft)

    print('nonzero k %s'%str(nonzero_k))
    pp(nonzero_v,msg='nonzero values')
    import amp_tools
    reload(amp_tools)
    values={}
    values['KA']=KA 
    values['KB']=KB
    values['KC']=KC
    values['A0']=A0
    values['A1']=A1
    values['B0']=B0
    values['B1']=B1
    values['C0']=C0
    values['C1']=C1
    values['Atheta']=Atheta
    values['Btheta']=Btheta
    values['Ctheta']=Ctheta
    expect = amp_tools.modes_1(**values)

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
    pp(expect_v,msg='expected values')

    ax0.scatter(expect_k,expect_v.real, label='xr',c='c')
    ax0.scatter(expect_k,expect_v.imag, label='xi',c='r')
    output.update({'Q_flat':Q_flat,'actual_fft':actual_fft,'qfft':qfft})
    output.update(values)
    output['expect']=expect
    output['nonzero_v']=nonzero_v
    output['nonzero_k']=nonzero_k
    output['expect_v']=expect_v
    output['expect_k']=expect_k
    output['B']=B
    if expect_k.size != nonzero_k[0].size:
        print("K sizes don't match")
    else:
        if np.abs( expect_k-nonzero_k).sum() > 1e-7:
            print("k values don't match!!!")
        print( "Total difference: %0.2e"%( np.abs( expect_v-nonzero_v).sum()))

    ax0.legend(loc=0)
    ax1.scatter(expect_k, (expect_v-nonzero_v).imag,label='imag',c='r')
    ax1.scatter(expect_k, (expect_v-nonzero_v).real,label='real',c='b')
    if np.abs(expect_k-nonzero_k).sum() > 0:
        print("VERY TERRIBLE ERROR K ERROR")
        ax1.set_title("SHIT K ERROR")
    ax1.set_yscale('symlog', linthreshy=1e-8)

    fig.savefig('p49c_plots/fft_1d_test2.png')
    plt.close(fig)
    return output
vals_small={   "A0" : 1, "A1" : np.sqrt(2.e-6), "KA":1, "Atheta" : 0,
   "B0" : 1, "B1" : np.sqrt(2.e-6), "KB":1, "Btheta" : 0, 
   "C0" : 1, "C1" : 0, "KC":5, "Ctheta":1.2}

vals_0={"A0" : 8, "A1" : 2, "KA":1, "Atheta" : 0,
        "B0" : 2, "B1" : 4, "KB":1, "Btheta" : 0, 
        "C0" : 1, "C1" : 0, "KC":5, "Ctheta":1.2}
vals_1={"A0" : 0, "A1" : 3, "KA":1, "Atheta" : 0,
        "B0" : 0, "B1" : 2, "KB":1, "Btheta" : 0, 
        "C0" : 1, "C1" : 0, "KC":5, "Ctheta":1.2}
oned=stuffs(**vals_1)
