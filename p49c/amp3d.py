if 0:
    #plt.legend(loc=0)
    #plt.savefig('fft7_test2.png')

    print('nonzero k %s'%str(nonzero_k))
    pp(nonzero_v,msg='nonzero values')
    def addit(exp,signed_k, label,v,phase_angle=0):
        k = tuple(np.abs(signed_k))
        phase_sign = np.sign(signed_k[-1])
        exp['labs'] = expect.get('labs',{})

        if np.abs(v) > 1e-10:
            exp[k] = exp.get(k,0)+v*np.exp(1j*phase_angle*phase_sign)
            exp['labs'][label]=k
            print('wut',v)#exp[k])

    expect={}

    addit( expect, np.zeros(3),       "0",                 A0*B0*C0)
    addit( expect, KA,      "KA",        C0*B0*(0.5*A1), Atheta)
    addit( expect, KB,      "KB",        A0*C0*(0.5*B1), Btheta)
                                      
    addit( expect, KA+KB,   "KA+KB",     C0*(0.5*A1)*(0.5*B1), Atheta+Btheta)
    addit( expect, KA-KB,   "KA-KB",     C0*(0.5*A1)*(0.5*B1), Atheta-Btheta)
                                      
"""
    addit( expect, KC,      "KC",        A0*B0*(0.5*C1), Ctheta)
    addit( expect, KB+KC,   "KB+KC",     A0*(0.5*B1)*(0.5*C1), Btheta+Ctheta)
    addit( expect, KB-KC,   "KB-KC",     A0*(0.5*B1)*(0.5*C1), Btheta-Ctheta)
                                      
    addit( expect, KA-KC,   "KA-KC",     B0*(0.5*A1)*(0.5*C1), (Atheta-Ctheta))
    addit( expect, KA+KC,   "KA+KC",     B0*(0.5*A1)*(0.5*C1), Atheta+Ctheta)
                                      
    addit( expect, KA+KB+KC,"KA+KB+KC",  (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta+Btheta+Ctheta)
    addit( expect, KA-KB+KC,"KA-KB+KC",  (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta-Btheta+Ctheta)
    addit( expect, -KA+KB+KC,"-KA+KB+KC", (0.5*A1)*(0.5*B1)*(0.5*C1), -Atheta+Btheta+Ctheta)
    addit( expect, KA+KB-KC,"KA+KB-KC",  (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta+Btheta-Ctheta)

    labs = expect.pop('labs')
    yv = np.logspace(-5,1,len(labs))
    for nmode,l in enumerate(labs):
        my_y=yv[nmode]
        plt.text( 10,my_y,l)
        plt.plot(  [labs[l],10], [my_y,my_y],c='k')
    expect_k = nar(sorted(list(expect.keys())))
    expect_v = nar([ expect[ key] for key in expect_k])
    plt.scatter(expect_k,expect_v.real, label='xr',c='c')
    plt.scatter(expect_k,expect_v.imag, label='xi',c='r')
    if expect_k.size != nonzero_k[0].size:
        print("K sizes don't match")
    else:
        if np.abs( expect_k-nonzero_k).sum() > 1e-7:
            print("k values don't match!!!")
        print( "Total difference: %0.2e"%( np.abs( expect_v-nonzero_v).sum()))
"""

from go import *
if '/home/dcollins/repos/p49c/p49_eigenmodes' not in sys.path:
    sys.path.append('/home/dcollins/repos/p49c/p49_eigenmodes')
if './p49c' not in sys.path:
    sys.path.append('./p49c')

import amp_tools
reload(amp_tools)


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
def pp(Q,msg=''):
    print('===')
    print("%s R %s"%(msg,str(Q.real)))
    print("%s I %s"%(msg,str(Q.imag)))
def lmax(inarr,dim):
    """plot the max of a 3d *inarr* along *dim*"""
    arr=np.zeros_like(inarr)
    ok =  inarr>1e-13
    arr[ok] = inarr[ok]
    out = np.max(np.abs(arr),axis=dim)
    return out

def stuff_3d(save=None):
    """in 3d!"""
    #x = np.arange(Q.size)
    output={}
    size=32
    twopi=np.pi*2
    x,y,z = np.mgrid[0:1:1./size,0:1:1./size,0:1:1./size]
    #x = np.mgrid[0:1:1./size]#,0:1:1./size,0:1:1./size]
    #y=z=1
    #y=0
    #z=0
    #k_unit_int = nar([1,7,1])
    #k_unit_Q = np.array(k_unit_int)*twopi
    #ampl=0.1
    #Q_flat =  ampl* (np.exp(1j*(k_unit_Q[0]*x+k_unit_Q[1]*y+k_unit_Q[2]*z))).real

    A0 = 1
    A1 = np.sqrt(2.e-6)
    KAx=0
    KAy=0
    KAz=1
    KA = nar([KAx,KAy,KAz])
    Aa = A0
    Atheta = 0
    Ab = A1* (np.exp(1j*(2*np.pi*(KAx*x+KAy*y+KAz*z)+Atheta))).real
    A =  Aa+Ab

    B0 = 1
    B1 = np.sqrt(2.e-6)
    KBx=0
    KBy=0
    KBz=1
    KB = nar([KBx,KBy,KBz])
    Ba = B0
    Btheta = 0 ;-np.pi/2.
    Bb = B1* (np.exp(1j*(2*np.pi*(KBx*x+KBy*y+KBz*z)+Btheta))).real
    B =  Ba+Bb
    output['B']=B

    C0 = 1
    C1 = 0
    KCx=5
    KCy=0
    KCz=0
    KC = nar([KCx,KCy,KCz])
    Ca = C0
    Ctheta=1.2
    Cb = C1* (np.exp(1j*(2*np.pi*(KCx*x+KCy*y+KCz*z)+Ctheta))).real
    C =  Ca+Cb

    plt.clf()
    slc = slice(3,4),slice(3,4),slice(None);the_x = z
    #slc = slice(None), slice(3,4),slice(3,4);the_x=x
    #slc = slice(None), slice(3,4),slice(3,4);the_x=x
    #slc = slice(None);the_x=x
    def xx(arr):
        return arr[slc].flatten()


    plt.plot(xx(the_x),xx(A),label='A')
    plt.plot(xx(the_x),xx(B),label='B')
    plt.plot(xx(the_x),xx(C),label='C')
    plt.plot(xx(the_x),xx(A*B*C),label='A*B*C')
    plt.legend(loc=0)
    plt.savefig('p49c_plots/fft7_test.png')
    if save is None:
        Q_flat = (A*B*C)
        actual_fft=np.fft.rfftn(Q_flat)
    else:
        Q_flat = (A*B*C)[slc]
        Q_flat = save['Q_flat']
        Q_flat.shape = Q_flat.size
        actual_fft=np.fft.rfft(Q_flat)
    qfft=actual_fft/(Q_flat.size)
    nonzero_k= nz(qfft)
    nonzero_v = nonzero(qfft)
    turb_quan.plotter2( [lmax( qfft.real, 0), 
                         lmax( qfft.imag, 0),
                         lmax( qfft.real, 1), 
                         lmax( qfft.imag, 1),
                         lmax( qfft.real, 2), 
                         lmax( qfft.imag, 2)],
                       "p49c_plots/qfft_3da.png", 
                       labs=['real x','imag x','real y','imag y','real z','imag z'],
                       share=False)
    def plotter3(arr,axis):
        dims = list(arr.shape)
        dims.pop(axis)
        flat = np.sum(arr,axis=axis)
        flat.shape=dims
        return flat
    turb_quan.plotter2( [plotter3( A, 0), 
                         plotter3( B, 0),
                         plotter3( C, 0), 
                         plotter3( A*B*C, 0)],
                       "p49c_plots/real_x_3da.png", 
                       labs=["Ax","Bx","Cx","ABCx"],
                       share=False,norm='ind')
    turb_quan.plotter2( [plotter3( A, 1), 
                         plotter3( B, 1),
                         plotter3( C, 1), 
                         plotter3( A*B*C, 1)],
                       "p49c_plots/real_y_3da.png", 
                       labs=["Ay","By","Cy","ABCy"],
                       share=False,norm='ind')
    turb_quan.plotter2( [plotter3( A, 2), 
                         plotter3( B, 2),
                         plotter3( C, 2), 
                         plotter3( A*B*C, 2)],
                       "p49c_plots/real_z_3da.png", 
                       labs=["Az","Bz","Cz","ABCz"],
                       share=False,norm='ind')


    values={}
    values['KA']=KA[2]
    values['KB']=KB[2]
    values['KC']=KC[2]
    values['A0']=A0
    values['A1']=A1
    values['B0']=B0
    values['B1']=B1
    values['C0']=C0
    values['C1']=C1
    values['Atheta']=Atheta
    values['Btheta']=Btheta
    values['Ctheta']=Ctheta
    reload(amp_tools)
    expect = amp_tools.modes_1(**values)
    labs = expect.pop('labs')
    output.update({'Q_flat':Q_flat,'acutal_fft':actual_fft,'qfft':qfft})
    output.update(values)
    expect_k = nar(sorted(list(expect.keys())))
    expect_v = nar([ expect[ key] for key in expect_k])
    output['expect']=expect
    output['nonzero_v']=nonzero_v
    output['nonzero_k']=nonzero_k
    output['expect_v']=expect_v
    output['expect_k']=expect_k
    return output
if 0:

    labs = expect.pop('labs')
    yv = np.logspace(-5,1,len(labs))
    for nmode,l in enumerate(labs):
        my_y=yv[nmode]
        plt.text( 10,my_y,l)
        plt.plot(  [labs[l],10], [my_y,my_y],c='k')
    pdb.set_trace()
    expect_k = nar(sorted(list(expect.keys())))
    expect_v = nar([ expect[ key] for key in expect_k])
    if expect_k.size != nonzero_k[0].size:
        print("K sizes don't match")
    else:
        if np.abs( expect_k-nonzero_k).sum() > 1e-7:
            print("k values don't match!!!")
        print( "Total difference: %0.2e"%( np.abs( expect_v-nonzero_v).sum()))

    plt.clf()
    #plt.plot(qfft.real,label='r')
    #plt.plot(qfft.imag,label='i')
    nonzero_k= nz(qfft)
    nonzero_v = nonzero(qfft)
    output={'Q_flat':Q_flat,'actual_fft':actual_fft,'qfft':qfft}
    output.update(values)
    output['expect']=expect
    output['nonzero_v']=nonzero_v
    output['nonzero_k']=nonzero_k
    output['expect_v']=expect_v
    output['expect_k']=expect_k
tred = stuff_3d()#save=oned)

def comp(a,b):
    for n in range(a.size):
        ar = a[n].real
        br = b[n].real
        print("%d %0.2e - %0.2e = %0.2e"%(n,ar,br,1-ar/br))

comp(tred['nonzero_v'],tred['expect_v'])
