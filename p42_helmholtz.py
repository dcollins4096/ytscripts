"""

Why were all the ghost zones set to -1?

"""
mark_time = None
import fourier_tools.fourier_filter as Filter
class CubeException(Exception):
    def __init__(self,filename):
        self.value = "Needs field %s"%filename
    def __str__(self):
        return repr(self.value)
    
def needs_fft(oober,frame, field_list,dtype='float32'):
    for field in field_list:
        filename = "%s/fft_%s.%s"%(oober.product_dir(frame),field,dtype)
        if glob.glob(filename) == []:
            raise CubeException(filename)
def MakeHelmholz(oober,frame,field,debug=1,dtype='float32'):
    #mark_time=time_marker()
    if field == 'velocity':
        fieldlist=['%s-velocity'%s for s in 'xyz']
    if field == 'acceleration':
        fieldlist=['%s-acceleration'%s for s in 'xyz']
    if field == 'Driving':
        fieldlist=['DrivingField%s'%s for s  in '123']
    needs_fft(oober,frame,fieldlist)
    
    vhat = oober.fft(frame,fieldlist[0],debug=debug)
    nx = vhat.shape
    kvec = np.ogrid[0:nx[0],0:nx[1],0:nx[2]]
    kdotv = vhat*kvec[0]
            
    for cmpt, cmpt_name in enumerate(fieldlist[1:]):
        vhat = oober.fft(frame,cmpt_name,debug=debug)
        kdotv += vhat*kvec[cmpt+1]

    del vhat
    NormK = kvec[0]**2+kvec[1]**2+kvec[2]**2
    NormK[0,0,0]=1.0
    if dtype == 'float32':
        fft_dtype = 'complex64'
    elif dtype == 'float64':
        fft_dtype = 'complex128'
    for cmpt, cmpt_name in enumerate(fieldlist):
        this_set = kvec[cmpt]*kdotv/(NormK)
        filename = '%s/fft_converging-%s.%s'%(oober.product_dir(frame),cmpt_name,dtype)
        fptr = h5py.File(filename,'w')
        fptr.create_dataset('converging-%s'%(cmpt_name),this_set.shape,data=this_set)
        fptr.close()
        print "Created",filename
        vhat = oober.fft(frame,cmpt_name,debug=debug)
        this_set = vhat - this_set
        filename = '%s/fft_solenoidal-%s.%s'%(oober.product_dir(frame),cmpt_name,dtype)
        fptr = h5py.File(filename,'w')
        fptr.create_dataset('solenoidal-%s'%(cmpt_name),this_set.shape,data=this_set)
        fptr.close()
        print "Created",filename

    #ConvergingPower(oobername,frame,field,debug,dtype)
    
def HelmholzPower(oober,frame,field,debug=1,dtype='float32'):
    #mark_time=time_marker()
    if field == 'velocity':
        fieldlist=['%s-velocity'%s for s in 'xyz']
    if field == 'acceleration':
        fieldlist=['%s-acceleration'%s for s in 'xyz']
    if field == 'Driving':
        fieldlist=['DrivingField%s'%s for s  in '123']
    #mark_time = time_marker()
    if dtype == 'float32':
        fft_dtype = 'complex64'
    elif dtype == 'float64':
        fft_dtype = 'complex128'
    else:
        fft_type=dtype
    power = 0
    for i,x in enumerate('xyz'):
        if mark_time is not None:
            mark_time('Start loop %s'%x)
        debug = 200
        Vhat = oober.fft(frame,'converging-%s-%s'%(x,field),num_ghost_zones=-1,debug=debug,dtype=dtype)
        power += (Vhat.conjugate()*Vhat)
        if mark_time is not None:
            mark_time('power addition')
    shell_average(power,oober,frame,'converging-%s'%field,debug,mark_time)
    power = 0
    for i,x in enumerate('xyz'):
        if mark_time is not None:
            mark_time('Start loop %s'%x)
        debug = 200
        Vhat = oober.fft(frame,'solenoidal-%s-%s'%(x,field),num_ghost_zones=0,debug=debug,dtype=dtype)
        power += (Vhat.conjugate()*Vhat)
        if mark_time is not None:
            mark_time('power addition')
    shell_average(power,oober,frame,'solenoidal-%s'%field,debug,mark_time)

def shell_average(power,oober,frame,field,debug=1,mark_time=None):
    ff = Filter.FourierFilter(power)
    if mark_time is not None:
        mark_time('made filter object')
    power_1d = np.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
    if mark_time is not None:
        mark_time('shell averages')
    filename = "%s/power_%s.h5"%(oober.product_dir(frame), field)
    if debug>0:
        print "spectra saved in ", filename
    file = h5py.File(filename,'w')
    file.create_dataset('power',power_1d.shape,data=power_1d)
    kspace=ff.get_shell_k()
    file.create_dataset('k',kspace.shape,data=kspace)
    if mark_time is not None:
        mark_time('saved power')
    file.close

def spectra_filename(oober,frame,xfield,field):
    dirname = oober.product_dir(frame)
    setname = 'power_%s.h5'%field
    if field in ['Density','LogDensity','gx','gy','gz']:
        setname = 'power_%s-work.h5'%field
    outname= "%s/%s"%(dirname,setname)
    print outname
    return outname

def MinK(TheY):
    TheY /= (TheY[ TheY != 0]).min()
    return TheY


def plot_helm(oober,frame,field):
    TheXmultiplier = MinK
    TheX_Lim=(9e-1,300 )
    TheX_Setname = 'k'
    TheY_Setname = 'power'
    TheX_Label=r'$k/k_{\rm{min}}$'
    fieldlist=['converging-%s'%field,'solenoidal-%s'%field]
    TheWeight = None
    plt.clf()
    out_power = []
    def do_log(f):
        return np.log10(f)
        #return f
    for n,field in enumerate(fieldlist):
        filename = spectra_filename(oober,frame,None,field)
        k,p = dpy(filename,['k','power'])
        out_power.append(p)
        plt.plot(do_log(k), do_log(p), label=['ud','us'][n])
    fname = '%s_%04d_Helmholtz_%s_ud_us.pdf'%(oober.outname,frame,field)
    plt.legend(loc=0)
    plt.savefig(fname)
    print fname
    return k,out_power[0], out_power[1]

