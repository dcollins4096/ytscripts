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

if 0:
    size=32
    twopi=np.pi*2
    x = np.mgrid[0:1:1./size]
    y=0
    z=0
    k_unit_int = nar([1,7,1])
    k_unit_Q = np.array(k_unit_int)*twopi
    ampl=0.1
    Q_flat =  ampl* (np.exp(1j*(k_unit_Q[0]*x+k_unit_Q[1]*y+k_unit_Q[2]*z))).real
    plt.clf()
    plt.plot(x,Q_flat)
    plt.savefig('fft_test.png')
    actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(0.5*x.size)

    pp(dumb(Q_flat,0),msg='k=0')
    pp(actual_fft[0],msg='k=0,actual')
    pp(dumb(Q_flat,1)/(actual_fft[1]),msg='k=1')
    pp(dumb(Q_flat,2),msg='k=2')
    pp(actual_fft[2],msg='k=2,actual')


if 0:
    plt.clf()
    plt.plot(qfft.real,label='r')
    plt.plot(qfft.imag,label='i')
    plt.scatter(k_unit_int[0], ampl,label='expect')
    plt.legend(loc=0)
    plt.savefig('fft_test2.png')

if 0:
    offset = 75
    Q_flat =  ampl* (np.exp(1j*(k_unit_Q[0]*x+k_unit_Q[1]*y+k_unit_Q[2]*z))).real
    Q_flat = Q_flat + offset
    plt.clf()
    plt.plot(x,Q_flat)
    plt.savefig('fft2_test.png')

    actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(x.size)

    plt.clf()
    plt.plot(qfft.real,label='r')
    nonzero_k= nz(qfft)
    nonzero_v = nonzero(qfft)
    expect_k = nar([0, k_unit_int[0]])
    expect_v = nar([offset, 0.5*ampl])
    print("error real %s"%str( (nonzero_v-expect_v).real))
    print("error imag %s"%str( (nonzero_v-expect_v).imag))

    plt.scatter( expect_k,expect_v,label='expect')
    plt.yscale('symlog', linthreshy=1e-5)
    plt.plot(qfft.imag,label='i')
    plt.legend(loc=0)
    plt.savefig('fft2_test2.png')

if 0:
    actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(x.size)
    ampl=7
    Q_1 =  ampl* (np.exp(1j*(k_unit_Q[0]*x+k_unit_Q[1]*y+k_unit_Q[2]*z))).real
    Q_flat = Q_1**2
    plt.clf()
    plt.plot(x,Q_flat)
    plt.savefig('fft3_test.png')
    actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(x.size)

    plt.clf()
    plt.plot(qfft.real,label='r')
    nonzero_k= nz(qfft)
    nonzero_v = nonzero(qfft)
    expect_k = nar([0, 2*k_unit_int[0]])
    expect_v = nar([0.5*ampl*ampl, (0.5*ampl)**2])

    for k in [0,1]:
        pp(dumb(Q_flat,k),msg='k=%d'%k)
        pp(actual_fft[k],msg='k=%d,actual'%k)

    print("expect %s"%str(expect_v))
    print("nonzero r %s"%str(nonzero_v.real))
    print("nonzero i %s"%str(nonzero_v.imag))
    print("error real %s"%str( (nonzero_v-expect_v).real))
    print("error imag %s"%str( (nonzero_v-expect_v).imag))

    plt.scatter( expect_k,expect_v,label='expect')
    plt.yscale('symlog', linthreshy=1e-5)
    plt.plot(qfft.imag,label='i')
    plt.legend(loc=0)
    plt.savefig('fft3_test2.png')

if 0:
    actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(x.size)
    ampl=7
    offset = 3
    Q_0 = offset
    Q_1 =  ampl* (np.exp(1j*(k_unit_Q[0]*x))).real
    Q_a = 2*Q_0*Q_1
    Q_b = Q_1**2
    Q_flat = (Q_0+Q_1)**2
    plt.clf()
    plt.plot(x,Q_flat)
    plt.plot(x,Q_0**2+Q_a+Q_b)
    plt.savefig('fft4_test.png')
    actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(x.size)

    plt.clf()
    plt.plot(qfft.real,label='r')
    nonzero_k= nz(qfft)
    nonzero_v = nonzero(qfft)
    print(nonzero_k)
    pp(nonzero_v,msg='ug')
    expect_k = nar([0, k_unit_int[0], 2*k_unit_int[0]])
    expect_v = nar([offset**2+0.5*ampl*ampl,0.5*2*offset*ampl, (0.5*ampl)**2])
    print("expect")
    print(expect_k)
    pp(expect_v,"exp")

if 0:
    for k in [0,1]:
        pp(dumb(Q_flat,k),msg='k=%d'%k)
        pp(actual_fft[k],msg='k=%d,actual'%k)

    print("expect %s"%str(expect_v))
    print("nonzero r %s"%str(nonzero_v.real))
    print("nonzero i %s"%str(nonzero_v.imag))
    print("error real %s"%str( (nonzero_v-expect_v).real))
    print("error imag %s"%str( (nonzero_v-expect_v).imag))

    plt.scatter( expect_k,expect_v,label='expect')
    plt.yscale('symlog', linthreshy=1e-5)
    plt.plot(qfft.imag,label='i')
    plt.legend(loc=0)
    plt.savefig('fft3_test2.png')

if 0:
    A0 = 3
    A1 = 4
    KA=1
    Aa = A0
    Ab = A1* (np.exp(1j*(2*np.pi*KA*x))).real
    A =  Aa+Ab
    B0 = 2
    B1 = 1.2
    KB=3
    Ba = B0
    Bb = B1* (np.exp(1j*(2*np.pi*KB*x))).real
    B =  Ba+Bb

    plt.clf()
    plt.plot(x,A,label='A')
    plt.plot(x,B,label='B')
    plt.plot(x,A+B,label='A+B')
    plt.legend(loc=0)
    plt.savefig('fft5_test.png')
    Q_flat = A*B
    actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(x.size)

    plt.clf()
    plt.plot(qfft.real,label='r')
    plt.plot(qfft.imag,label='i')
    plt.legend(loc=0)
    plt.yscale('symlog', linthreshy=1e-5)
    expect_k=nar([0,    KA,np.abs(KA-KB),KB,KA+KB])
    args = np.argsort(expect_k)

    expect_v=nar([A0*B0, 0.5*B0*A1, 0.25*A1*B1, 0.5*B1*A0, 0.25*A1*B1])
    plt.scatter(expect_k[args],expect_v[args],label='x')
    nonzero_k= nz(qfft)
    nonzero_v = nonzero(qfft)
    plt.savefig('fft5_test2.png')
    print('nonzero k %s'%str(nonzero_k))
    pp(nonzero_v,msg='nonzero values')
    pp(expect_v[args]-nonzero_v,msg='diff')
#   expect_k = nar([0, k_unit_int[0], 2*k_unit_int[0]])
#   expect_v = nar([offset**2+0.5*ampl*ampl,0.5*2*offset*ampl, (0.5*ampl)**2])
#   print("expect")
#   print(expect_k)
#   pp(expect_v,"exp")
if 0:
    A0 = 0
    A1 = 4
    KA=1
    Aa = A0
    Ab = A1* (np.exp(1j*(2*np.pi*KA*x))).real
    A =  Aa+Ab
    B0 = 0
    B1 = 1.2
    KB=3
    Ba = B0
    Bb = B1* (np.exp(1j*(2*np.pi*KB*x))).real
    B =  Ba+Bb

    C0 = 0
    C1 = 4.5
    KC=5
    Ca = C0
    Cb = C1* (np.exp(1j*(2*np.pi*KC*x))).real
    C =  Ca+Cb

    plt.clf()
    plt.plot(x,A,label='A')
    plt.plot(x,B,label='B')
    plt.plot(x,C,label='C')
    plt.plot(x,A*B*C,label='A*B*C')
    plt.legend(loc=0)
    plt.savefig('fft6_test.png')
    Q_flat = A*B*C
    actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(x.size)

    plt.clf()
    plt.plot(qfft.real,label='r')
    plt.plot(qfft.imag,label='i')
    plt.legend(loc=0)
    plt.yscale('symlog', linthreshy=1e-5)
    nonzero_k= nz(qfft)
    nonzero_v = nonzero(qfft)

    print('nonzero k %s'%str(nonzero_k))
    pp(nonzero_v,msg='nonzero values')

    expect_k=nar([np.abs(KA+KB+KC),np.abs(KA+KB-KC),np.abs(KA-KB+KC),np.abs(KB-KA+KC)])
    print(sorted(expect_k))
#   args = np.argsort(expect_k)

#   expect_v=nar([A0*B0, 0.5*B0*A1, 0.25*A1*B1, 0.5*B1*A0, 0.25*A1*B1])
#   plt.scatter(expect_k[args],expect_v[args],label='x')

#   pp(expect_v[args]-nonzero_v,msg='diff')
    plt.savefig('fft6_test2.png')


if 1:

    A0 = 1
    A1 = 4
    KA=1
    Aa = A0
    Ab = A1* (np.exp(1j*(2*np.pi*KA*x))).real
    A =  Aa+Ab

    B0 = 8
    B1 = 1.2
    KB=3
    Ba = B0
    Bb = B1* (np.exp(1j*(2*np.pi*KB*x))).real
    B =  Ba+Bb

    C0 = 0
    C1 = 4.5
    KC=5
    Ca = C0
    Cb = C1* (np.exp(1j*(2*np.pi*KC*x))).real
    C =  Ca+Cb

    plt.clf()
    plt.plot(x,A,label='A')
    plt.plot(x,B,label='B')
    plt.plot(x,C,label='C')
    plt.plot(x,A*B*C,label='A*B*C')
    plt.legend(loc=0)
    plt.savefig('fft6_test.png')
    Q_flat = A*B*C
    actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(x.size)

    plt.clf()
    plt.plot(qfft.real,label='r')
    plt.plot(qfft.imag,label='i')
    plt.legend(loc=0)
    plt.yscale('symlog', linthreshy=1e-5)
    nonzero_k= nz(qfft)
    nonzero_v = nonzero(qfft)

    print('nonzero k %s'%str(nonzero_k))
    pp(nonzero_v,msg='nonzero values')
    def addit(exp,k,v):
        exp[k] = exp.get(k,0)+v

    expect={}
    addit( expect, np.abs(KA+KB+KC), (0.5*A1)*(0.5*B1)*(0.5*C1))
    addit( expect, np.abs(KA-KB+KC), (0.5*A1)*(0.5*B1)*(0.5*C1))
    addit( expect, np.abs(-KA+KB+KC), (0.5*A1)*(0.5*B1)*(0.5*C1))
    addit( expect, np.abs(KA+KB-KC), (0.5*A1)*(0.5*B1)*(0.5*C1))

    addit( expect, np.abs(KB+KC), A0*(0.5*B1)*(0.5*C1))
    addit( expect, np.abs(KB-KC), A0*(0.5*B1)*(0.5*C1))

    addit( expect, np.abs(KA-KC), B0*(0.5*A1)*(0.5*C1))
    addit( expect, np.abs(KA+KC), B0*(0.5*A1)*(0.5*C1))
    addit( expect, np.abs(KC), A0*B0*(0.5*C1))
    #expect[ 0+0+0] = A0*B0*C0
    #expect[KA+KB+KC] = (0.5*A1)*(0.5*B1)*(0.5*C1)
    #expect[np.abs(KA+KB-KC)] = (0.5*A1)*(0.5*B1)*(0.5*C1)
    expect_k = nar(sorted(list(expect.keys())))
    expect_v = nar([ expect[ key] for key in expect_k])
    plt.scatter(expect_k,expect_v, label='x')

    #expect_k=nar([np.abs(KA+KB+KC),np.abs(KA+KB-KC),np.abs(KA-KB+KC),np.abs(KB-KA+KC),
    #             np.abs(KB),np.abs(KC),np.abs(KB+KC),np.abs(KB-KC)])
    #print(sorted(expect_k))
#   args = np.argsort(expect_k)

#   expect_v=nar([A0*B0, 0.5*B0*A1, 0.25*A1*B1, 0.5*B1*A0, 0.25*A1*B1])
#   plt.scatter(expect_k[args],expect_v[args],label='x')

#   pp(expect_v[args]-nonzero_v,msg='diff')
    plt.savefig('fft6_test2.png')
