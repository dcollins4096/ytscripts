execfile('go')
frame = 60; basedir = '/scratch1/dcollins/Paper08/B02/512/'
fname = '%s/RS%04d/restart%04d'%(basedir, frame,frame)
peaks = fPickle.load('p12_b02_512_n0060_peaks2.pickle')

ds = yt.load(fname)
if 'mhd' not in dir():
    ef('p13_fields.py')
fields = ['Lambda_therm_pos', 'Lambda_therm_neg', 'Lambda_turb_pos','Lambda_turb_neg',
          'Lambda_magn_pos', 'Lambda_magn_neg']
Colors = {'Lambda_therm_pos':'b', 'Lambda_therm_neg':'b', 'Lambda_turb_pos':'r','Lambda_turb_neg':'r',
          'Lambda_magn_pos':'g', 'Lambda_magn_neg':'g'}
Lines  = {'Lambda_therm_pos':'-', 'Lambda_therm_neg':'--', 'Lambda_turb_pos':'-','Lambda_turb_neg':'--',
          'Lambda_magn_pos':'-', 'Lambda_magn_neg':'--'}
Gravity = []
AllPlots = {}
for field in fields:
    AllPlots[field]=[]

all_fields = ["CellVolumeNorm", "Delta", "Gravity", "KineticEnergyDensity", "Pressure", 
       "GradProductDensityPressure", "LaplacePressure", "Lambda_therm", "Lambda_therm_pos", 
       "Lambda_therm_neg", "EnstrophyDensity", "RateOfStrainSqr", "Denstrophy", "DiffOmegaSqrStrainSqr", 
       "Lambda_turb", "Lambda_turb_pos", "Lambda_turb_neg", "MagneticPressure", "GradProductDensityMagneticPressure", 
       "LaplaceMagneticPressure", "MagneticFieldCrossContractions", "Lambda_magn", "Lambda_magn_pos", "Lambda_magn_neg"]

t0 = time.time()

if 1:
    """for checking the bad oklahoma"""
    peak = peaks[0]
    width = 0.05
    data = ds.sphere(peak,(width,'code_length'))
    field = 'velocity_divergence'
    #field = 'Lambda_magn_neg'
    ax=0
    proj = ds.proj(field, ax, data_source = data, center = peak)
    pw = proj.to_pw(center = peak, width = 2*width)
    pw.set_zlim('velocity_divergence',-100,100)
    pw.set_log('velocity_divergence',False)

    pw.save('yt3-glitch-smoothed-false')
    pw.save('yt3-glitch-mjt-fix')

if 0:
    fields = ['Lambda_therm','Lambda_turb','Lambda_magn','density','velocity_divergence']
    Fcol = {'Lambda_therm':'b','Lambda_turb':'r','Lambda_magn':'g','density':'k','velocity_divergence':'c'}
    Fcol.update({'Lambda_therm_pos_per_delta':'b', 'Lambda_turb_pos_per_delta':'r','Lambda_turb_neg_per_delta':'r',
                 'Lambda_therm_neg_per_delta':'b','Lambda_magn_pos_per_delta':'g','Lambda_magn_neg_per_delta':'g',
                 'Lambda_therm_per_delta':'b','Lambda_magn_per_delta':'g','Lambda_turb_per_delta':'r'})
    Line = {'Lambda_turb_neg_per_delta':'--', 'Lambda_therm_neg_per_delta':'--','Lambda_magn_neg_per_delta':'--',
            'Lambda_turb_pos_per_delta':'-.', 'Lambda_therm_pos_per_delta':'-.','Lambda_magn_pos_per_delta':'-.'} 


    symlog_thresh = {'Lambda_therm':1,'Lambda_turb':1,'Lambda_magn':1,'velocity_divergence':1}
    labels = dict(zip(['Lambda_therm','Lambda_turb','Lambda_magn','density','velocity_divergence'],
                      ['LH','LT','LM','rho','d']))
    Fcol['Lambda_therm_per_delta']='b'
    labels['Lambda_therm_per_delta'] = 'LH/d'
    labels.update({'Lambda_therm_pos_per_delta':'LHp/d', 'Lambda_turb_pos_per_delta':'LTp/d','Lambda_turb_neg_per_delta':'LTn/d',
                 'Lambda_therm_neg_per_delta':'LHn/d','Lambda_magn_pos_per_delta':'LMp/d','Lambda_magn_neg_per_delta':'LMp/d',
                 'Lambda_therm_per_delta':'LHt/d','Lambda_magn_per_delta':'LMt/d','Lambda_turb_per_delta':'LTt/d'})
    Rlist = [0.5e-3, 1e-3, 1e-2, 5e-2]
    Rcolo = [   'r',  'k',  'r',  'k']

    if 1:
        these_fields=['Lambda_therm_pos_per_delta', 'Lambda_turb_pos_per_delta','Lambda_turb_neg_per_delta',
                     'Lambda_therm_neg_per_delta','Lambda_magn_pos_per_delta','Lambda_magn_neg_per_delta',
                      'Lambda_therm_per_delta', 'Lambda_turb_per_delta', 'Lambda_magn_per_delta'][0:1]
        if 0:
            data = ds.all_data()
        if 1:
            """for testing"""
            width = 0.05
            peak = peaks[0]
            data = ds.sphere(peak,(width,'code_length'))

        print "adding", these_fields[0:1]
        prof = yt.create_profile(data,'density',fields=these_fields[0:1], weight_field='cell_volume')
        plt.clf()
        for field in these_fields[1:]:
            print "adding", field
            prof.add_fields(field)
        for field in these_fields:
            plt.plot(0.5*(prof.x_bins[1:]+prof.x_bins[:-1]), prof[field],linestyle = Line.get(field,'-'), label=labels[field], c=Fcol[field])
        plt.yscale('symlog')
        plt.xscale('log')
        plt.legend(loc=0)
        fname = 'b025v_full'
        print fname
        plt.savefig(fname)

    if 0:
        """
        adding Lambda_turb_pos_per_delta
        adding Lambda_turb_neg_per_delta
        adding Lambda_therm_neg_per_delta
        adding Lambda_magn_pos_per_delta
        adding Lambda_magn_neg_per_delta
        adding Lambda_therm_per_delta
        adding Lambda_turb_per_delta
        adding Lambda_magn_per_delta
        #for fields in these_fields[1:]:
            """


    if 0:
        # for science.
        for npeak, peak in enumerate(peaks):
            print "peak", npeak
            width = 0.05
            data = ds.sphere(peak,(width,'code_length'))
            ax=0

            if 0:
                for field in fields:
                    proj = ds.proj(field, ax, data_source = data, center = peak)
                    pw = proj.to_pw(center = peak, width = 2*width)
                    for ncol,R in enumerate(Rlist):
                        pw.annotate_sphere(peak, (R,'code_length'), circle_args={'color':Rcolo[ncol]})
                    if symlog_thresh.has_key(field):
                        pw.set_log(field,True,linthresh=symlog_thresh[field])
                    pw.save('b025_n%04d_peak%04d'%(frame,npeak))
            if 0:
                prof=yt.create_profile(data, 'radius',fields=fields,weight_field='cell_volume')

            plt.clf()
            for field in fields:
                plt.plot(0.5*(prof.x_bins[1:]+prof.x_bins[:-1]), prof[field],label=labels[field], c=Fcol[field])

            ylim = plt.ylim()
            for ncol,R in enumerate(Rlist):
                plt.plot([R,R],ylim,linestyle='dotted', c=Rcolo[ncol]) 
            plt.ylim(ylim)
            plt.yscale('symlog')
            plt.legend(loc=0)

            plt.xscale('linear')
            fname = 'b025_n%04d_peak%04d_profiles_linear'%(frame,npeak)
            print fname
            plt.savefig(fname)
            plt.xscale('log')
            fname = 'b025_n%04d_peak%04d_profiles_log'%(frame,npeak)
            print fname
            plt.savefig(fname)


#plot = yt.ProfilePlot.from_profiles(prof)
#print plot.save()
            
