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

def stuff_3d(A=[],B=[],K=[],Theta=[],linthreshy=1e-8):
    """in 3d!"""
    #x = np.arange(Q.size)
    output={}
    size=32
    twopi=np.pi*2
    x,y,z = np.mgrid[0:1:1./size,0:1:1./size,0:1:1./size]

    nmodes = len(A)
    consts=np.zeros([nmodes])
    linears=np.zeros([nmodes,size,size,size])
    sums=np.zeros([nmodes,size,size,size])
    for n in range(nmodes):
        consts[n]=A[n]
        linears[n,...]=B[n]*(np.exp(1j*(2*np.pi*(K[n][0]*x+K[n][1]*y+K[n][2]*z)+Theta[n])))
        sums[n,...] = consts[n]+linears[n,...]

    plt.clf()
    slc = slice(3,4),slice(3,4),slice(None);the_x = z
    #slc = slice(None), slice(3,4),slice(3,4);the_x=x
    #slc = slice(None), slice(3,4),slice(3,4);the_x=x
    #slc = slice(None);the_x=x
    def xx(arr):
        return arr[slc].flatten()


    for n in range(nmodes):
        plt.plot(xx(the_x),xx(sums[n,...]),label=n)
    Q_flat = np.prod(sums.real,axis=0)
    plt.plot(xx(the_x),xx(Q_flat),label='A*B*C')
    plt.legend(loc=0)
    plt.savefig('p49c_plots/fft_3db_test.png')
    actual_fft=np.fft.rfftn(Q_flat)
    qfft=actual_fft/(Q_flat.size)
    nonzero_k= nz(qfft)
    nonzero_v = nonzero(qfft)
    abs_qfft =np.abs(qfft)
    zmin = abs_qfft[ abs_qfft > 0].min()
    turb_quan.plotter2( [lmax( qfft.real, 0), 
                         lmax( qfft.imag, 0),
                         lmax( qfft.real, 1), 
                         lmax( qfft.imag, 1),
                         lmax( qfft.real, 2), 
                         lmax( qfft.imag, 2)],
                       "p49c_plots/qfft_3da.png", 
                       norm='positive',zmin=zmin,
                       labs=['real x','imag x','real y','imag y','real z','imag z'],
                       share=False)
    def plotter3(arr,axis):
        dims = list(arr.shape)
        dims.pop(axis)
        flat = np.sum(arr,axis=axis)
        flat.shape=dims
        return flat
    turb_quan.plotter2( [plotter3( sums[0,...], 0), 
                         plotter3( sums[1,...], 0),
                         plotter3( sums[0,...], 0), 
                         plotter3( Q_flat, 0)],
                       "p49c_plots/real_x_3da.png", 
                       labs=["Ax","Bx","Cx","ABCx"],
                       share=False,norm='ind')
    turb_quan.plotter2( [plotter3( sums[0,...], 1), 
                         plotter3( sums[1,...], 1),
                         plotter3( sums[0,...], 1), 
                         plotter3( Q_flat, 1)],
                       "p49c_plots/real_y_3da.png", 
                       labs=["Ay","By","Cy","ABCy"],
                       share=False,norm='ind')
    turb_quan.plotter2( [plotter3( sums[0,...], 2), 
                         plotter3( sums[1,...], 2),
                         plotter3( sums[0,...], 2), 
                         plotter3( Q_flat, 2)],
                       "p49c_plots/real_z_3da.png", 
                       labs=["Az","Bz","Cz","ABCz"],
                       share=False,norm='ind')


    values={}
    values['A'] = nar(A)
    values['B'] = nar(B)
    values['K'] = nar(K)
    values['Theta']=nar(Theta)
    reload(amp_tools)
    expect = amp_tools.modes_2(**values)
    output['labs']= expect.pop('labs')
    output['full']= expect.pop('full_history')

    output.update({'Q_flat':Q_flat,'acutal_fft':actual_fft,'qfft':qfft})
    output.update(values)
    output['expect']=expect


    output['nonzero_v']=nonzero_v
    output['nonzero_k']=zip(*nonzero_k)
    output['nonz']= dict( zip(output['nonzero_k'], output['nonzero_v']))
    def sanitize_real_vectors(ft):
        for v in ft:
            if v[2] == 0:
                pass
    #sanitize_real_vectors( output['nonz'])
    total_error = 0
    for n,vec in enumerate(output['nonzero_k']):
        if vec in expect:
            this_error = np.abs( expect[vec] - nonzero_v[n])
            total_error += this_error
            if this_error > 1e-12:
                print("expect off. actual %0.2e expect %0.2e diff %0.2d vec %s"%\
                      ( nonzero_v[n], expect[vec], this_error, str(vec)))
        else:
            print("expect missing vector %s"%str(vec))
    print("Total error: %0.2e"%total_error)

    #expect_k = nar(sorted(list(expect.keys())))
    #expect_v = nar([ expect[ key] for key in expect_k])
    #output['expect_v']=expect_v
    #output['expect_k']=expect_k
    output['total_error']=total_error
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
def comp(a,b):
    for n in range(a.size):
        ar = a[n].real
        br = b[n].real
        print("%d %0.2e - %0.2e = %0.2e"%(n,ar,br,1-ar/br))
xhat=nar([1,0,0])
yhat=nar([0,1,0])
zhat=nar([0,0,1])
#vals_n = {'A':nar([0,0]),'B':nar([1e-2,1e-2]),'K':nar([[0,0,1],[0,0,1]]),'Theta':nar([0,0])}
vals_n = {'A':nar([0,0,1]),'B':nar([1e-2,1e-2,0]),'K':nar([[0,0,1],[0,0,1],[0,0,0]]),'Theta':nar([0,0,0])}
#vals_n = {'A':nar([1,1]),'B':nar([1,1]),'K':nar([yhat+zhat,yhat+2*zhat]),'Theta':nar([0.,0.])}
#vals_n = {'A':nar([1,1]),'B':nar([1,0]),'K':nar([zhat,zhat]),'Theta':nar([0.0,0.0])}#nar([0,0.2,0.1,np.pi/2,0])}

vals_n = {'A':nar([1,1,1,1,1]),'B':nar([1,1,1,1,1]),'K':nar([zhat,2*zhat,yhat+zhat,zhat,zhat]),'Theta':nar([0.2,0.1,np.pi/2,0,1.3,  4])}#doesnt work yet
vals_n = {'A':nar([1,1]),'B':nar([1,1]),'K':nar([zhat,yhat+zhat]),'Theta':nar([0.2,0.1,np.pi/2,0,1.3,  4])}#doesnt work yet
vals_n = {'A':nar([1,1]),'B':nar([1,1]),'K':nar([zhat+xhat,yhat+2*zhat]),'Theta':nar([0.,0.])}#doesnt work yet
stre=stuff_3d(**vals_n)
