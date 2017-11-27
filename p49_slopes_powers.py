import matplotlib
import pdb
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import fPickle
import pylab
from scipy.optimize import leastsq
import sys
import glob
plt.close('all')

def linear_chi2(par,x,y) :
    m = par[0]
    y0 = par[1]

    xcent = (x[0]+x[-1])/2

    ymodel = m*(x-xcent)+y0
#    ymodel = m*(x)+y0
    
    return(y - ymodel)




def slopes_powers(car,frame, n0=1,p=1, fitlmin=1e3, fitlmax=8e3, plot_format="pdf"):
    base = '%s/FRBs/DD%04d_Cl%s_n0-%04d_p-%d.dat'%(car.directory,frame, "%s", n0, p)
    output_base = '%s_DD%04d_n0-%04d_p-%d_%s.%s'%(car.outname,frame, n0, p,"%s","%s")
    shortname = car.name
    shortname_n0p = "%s_n0-%04d_p-%d"%(shortname, n0, p)

#fitlmin = 5e3
#fitlmax = 3e4


    ClE = {}
    ClB = {}

    for ax in ['x','y','z'] :
        ell, ClE[ax], ClB[ax] = pylab.loadtxt(base % ax,usecols=[0,1,2],unpack=True)
        
    ellmask = (fitlmin < ell)*(ell < fitlmax)

    Eslope = {}
    Bslope = {}

    Eamp = {}
    Bamp = {}

    import pdb
    #pdb.set_trace()
    for ax in ['x','y','z'] :
        mb = np.array([0,0])
        ellmask_finite = np.logical_and(ellmask, ClE[ax] > 0)
        ellmask_finite = np.logical_and(ellmask_finite, ClB[ax] > 0)


        if ellmask_finite.sum() == 0:
            make_plots=False
            print "Warning: No positive values of harmonics in target ell range"
            Bslope[ax]=np.nan
            Bamp[ax]=np.nan
            Eslope[ax]=np.nan
            Eamp[ax]=np.nan
            continue
        make_plots=True
        x =np.log10(ell[ellmask_finite])

        y = np.log10(ClE[ax][ellmask_finite])
        res = leastsq(linear_chi2, mb, args=(x,y) )
        Eslope[ax] = res[0][0]
        Eamp[ax] = pow(10,res[0][1])
        
        y = np.log10(ClB[ax][ellmask_finite])
        res = leastsq(linear_chi2, mb, args=(x,y) )
        Bslope[ax] = res[0][0]
        Bamp[ax] = pow(10,res[0][1])

        
    if make_plots:
        coldict = {'x':'r', 'y':'g', 'z':'b'}
        #values_dict = {}
        ##pickle_fname = 'p49_all_values.pickle'
        ##if glob.glob(pickle_fname) != []:
        ##    values_dict = fPickle.load(pickle_fname)
        #    #[shortname][frame][fit range, Eslope(3) Bslope(3) Eamp(3) Bamp(3)]
        #if not values_dict.has_key(shortname_n0p):
        #    values_dict[shortname_n0p]={}
        #if not values_dict[shortname_n0p].has_key(frame):
        #    values_dict[shortname_n0p][frame]={'fit_range':[fitlmin, fitlmax]}
        #    values_dict[shortname_n0p][frame]['Eslope'] = Eslope
        #    values_dict[shortname_n0p][frame]['Bslope'] = Bslope
        #    values_dict[shortname_n0p][frame]['Eamp'] = Eamp
        #    values_dict[shortname_n0p][frame]['Bamp'] = Bamp
        ##fPickle.dump(values_dict,pickle_fname)
            

# plots
        lfit = pow(10,(x[0]+x[-1])/2)

        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.plot(ell,Eamp['x']*(ell/lfit)**(-2.43), c=[0.5,0.5,0.5], label='%0.2f'%(-2.43))
        for ax in ['x','y','z'] :
            axis.plot(ell,ClE[ax],'-',label='E%s %0.2f'%(ax, Eslope[ax]), c=coldict[ax])
            axis.plot(ell,ClB[ax],'--',label='B%s %0.2f'%(ax, Bslope[ax]), c=coldict[ax])
            axis.plot(ell,Eamp[ax]*(ell/lfit)**Eslope[ax],'k:')
            axis.plot(ell,Bamp[ax]*(ell/lfit)**Bslope[ax],'k:')

            print('%s Eslope %f Bslope %f Eamp/Bamp %f' % (ax,Eslope[ax],Bslope[ax],Eamp[ax]/Bamp[ax]))
            
        axis.legend()
        axis.loglog()
        axis.set_xlabel(r'$\ell$')
        axis.axvspan(fitlmax,pylab.gca().get_xlim()[1],facecolor='k',alpha=0.1,ls=None)
        axis.axvspan(pylab.gca().get_xlim()[0],fitlmin,facecolor='k',alpha=0.1,ls=None)
        fig.savefig(output_base%("EE_BB",plot_format))

        fig = plt.figure()
        axis=fig.add_subplot(111)
        axis.plot(ell, np.zeros_like(ell)+2.0, c=[0.5]*3)
        for ax in ['x','y','z'] :
            axis.plot(ell,ClE[ax]/ClB[ax],'-',label=ax, c=coldict[ax])

        axis.set_title('ClE/ClB')
        axis.legend()
        axis.loglog()
        axis.set_xlabel(r'$\ell$')
        fig.savefig(output_base%("Ratio",plot_format))



        fig = plt.figure()
        axis=fig.add_subplot(111)
        ax = 'y'
        axis.plot(ell,ClE[ax],'-',label='E%s'%ax, c=coldict[ax])
        axis.plot(ell,ClB[ax],'--',label='B%s'%ax, c=coldict[ax])
        axis.plot(ell,Eamp[ax]*(ell/lfit)**Eslope[ax],'k:')
        axis.plot(ell,Bamp[ax]*(ell/lfit)**Bslope[ax],'k:')

        axis.set_xlim(500,1.3e4)

#legend()
        axis.loglog()

        axis.axvspan(pylab.gca().get_xlim()[0],fitlmin,facecolor='k',alpha=0.1,ls=None)
        axis.axvspan(fitlmax,pylab.gca().get_xlim()[1],facecolor='k',alpha=0.1,ls=None)

#fig.savefig('EBfit.pdf',bbox_inches='tight')
        axis.set_xlabel(r'$\ell$')
        fig.savefig(output_base%("EBfit",plot_format))

    return {'Eslope':Eslope,'Eamp':Eamp,'Bslope':Bslope,'Bamp':Bamp}
