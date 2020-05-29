
from go import *
def rb(kall, power, nbins = 32,bins=None):    
    kmag = np.sqrt(kall[0,...]**2+kall[1,...]**2+kall[2,...]**2)
    dbin = (kmag.max()-kmag.min())/nbins
    if bins is None:
        bins = np.linspace(kmag.min()-0.5*dbin,kmag.max()+0.5*dbin,nbins+1)
    else:
        nbins = bins.size-1
    bin_center = 0.5*(bins[1:]+bins[:-1])
    dig=np.digitize(kmag,bins)
    
    power_1d = np.array([ power[ dig==this_bin].sum() for this_bin in range(1,nbins+1)])
    k2 = np.array([ (dig==this_bin).sum() for this_bin in range(1,nbins+1)])
    power_1d/=k2
    pdb.set_trace()


    fig,ax=plt.subplots(1,1,figsize=(12,12))

    ax.plot(bin_center,power_1d,c='k')
    fig.savefig('p56_binned.png')




    return bin_center, power_1d

