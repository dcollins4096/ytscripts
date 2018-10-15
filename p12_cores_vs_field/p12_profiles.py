execfile('go')
frame = 60; basedir = '/scratch1/dcollins/Paper08/B02/512/'
fname = '%s/RS%04d/restart%04d'%(basedir, frame,frame)
peaks = fPickle.load('../yt3_scripts/p12_b02_512_n0060_peaks2.pickle')

ds = yt.load(fname)
#ef('Paper13/p13_fields_for_paper12.py')
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


fields = ['Lambda_therm','Lambda_turb','Lambda_magn','density','velocity_divergence']
Fcol = {'Lambda_therm':'b','Lambda_turb':'r','Lambda_magn':'g','density':'m','Density':'m','velocity_divergence':'c','Lambda':'k'}
Fcol.update({'Lambda_therm_pos_per_delta':'b', 'Lambda_turb_pos_per_delta':'r','Lambda_turb_neg_per_delta':'r',
             'Lambda_therm_neg_per_delta':'b','Lambda_magn_pos_per_delta':'g','Lambda_magn_neg_per_delta':'g',
             'Lambda_therm_per_delta':'b','Lambda_magn_per_delta':'g','Lambda_turb_per_delta':'r'})
Line = {'Lambda_turb_neg_per_delta':'--', 'Lambda_therm_neg_per_delta':'--','Lambda_magn_neg_per_delta':'--',
        'Lambda_turb_pos_per_delta':'-.', 'Lambda_therm_pos_per_delta':'-.','Lambda_magn_pos_per_delta':'-.'} 
Fcol['Gravity']=[1,105./256,0,1]


symlog_thresh = {'Lambda_therm':1,'Lambda_turb':1,'Lambda_magn':1,'velocity_divergence':1}
labels = dict(zip(['Lambda_therm','Lambda_turb','Lambda_magn','density','velocity_divergence'],
                  ['LH','LT','LM','rho','d']))
Fcol['Lambda_therm_per_delta']='b'
Fcol['Lambda_per_delta']='k'
labels['Lambda_therm_per_delta'] = 'LH/d'
labels['Density']='d'
labels['Gravity']='G'
labels['Lambda_per_delta']='L/d'
labels.update({'Lambda_therm_pos_per_delta':'LHp/d', 'Lambda_turb_pos_per_delta':'LTp/d','Lambda_turb_neg_per_delta':'LTn/d',
             'Lambda_therm_neg_per_delta':'LHn/d','Lambda_magn_pos_per_delta':'LMp/d','Lambda_magn_neg_per_delta':'LMp/d',
             'Lambda_therm_per_delta':'LHt/d','Lambda_magn_per_delta':'LMt/d','Lambda_turb_per_delta':'LTt/d',
               'Lambda':'L','Lambda_per_delta':'L/d'})
Rlist = [0.5e-3, 1e-3, 1e-2, 5e-2]
Rcolo = [   'r',  'k',  'r',  'k']

if 1:
    for npeak, peak in enumerate(peaks):
        if npeak not in [0]:
            continue
        these_fields=['Lambda_therm_pos_per_delta', 'Lambda_turb_pos_per_delta','Lambda_turb_neg_per_delta',
                     'Lambda_therm_neg_per_delta','Lambda_magn_pos_per_delta','Lambda_magn_neg_per_delta',
                      'Lambda_therm_per_delta', 'Lambda_turb_per_delta', 'Lambda_magn_per_delta']
        if 0:
            these_fields=['Lambda_therm','Lambda_magn','Lambda_turb','Lambda','Gravity']
            vert_label = r'$\Lambda$'
            file_label = 'Lambda'
        if 1:
            these_fields=['Lambda_therm_per_delta','Lambda_magn_per_delta','Lambda_turb_per_delta','Lambda_per_delta','Gravity']
            vert_label = r'$\Lambda/4 \pi G \rho_0 \delta$'
            file_label = 'Lambda_G'
        if 1:
            """for testing"""
            width = 0.05
            data = ds.h.sphere(peak,width)

        clobber=False
        for bin_field in ['density','radius']:
            fname = 'b025v_n%04d_peak%04d_%s_%s'%(frame,npeak,bin_field,file_label)
            picklename = 'b025v_n%04d_peak%04d_%s_%s.pickle'%(frame,npeak,bin_field,file_label)
            plt.clf()
            if glob.glob(picklename) != [] and clobber==False:
                print "reading stuff from ", picklename
                field_data = fPickle.load(picklename)
            else:
                print "adding", these_fields[0:1]
                prof = yt.create_profile(data,bin_field,fields=these_fields[0:1], weight_field='cell_volume')
                for field in these_fields[1:]:
                    print "adding", field
                    prof.add_fields(field)
                prof.field_data['x_bin']=prof.x_bins
                field_data = prof.field_data
                fPickle.dump(prof.field_data,picklename)
            for field in these_fields:
                #fff = prof.field_data[('gas',field)]
                fff = field_data[field]
                plt.plot(0.5*(field_data['x_bin'][1:]+field_data['x_bin'][:-1]), fff,linestyle = Line.get(field,'-'), label=labels[field], c=Fcol[field])
            if bin_field == 'Radius':
                """for R"""
                ylim = plt.ylim()
                for ncol,R in enumerate(Rlist):
                    plt.plot([R,R],ylim,linestyle='dotted', c=Rcolo[ncol])
                plt.ylim(ylim)
                plt.xscale('linear')
                plt.xlabel('R')
                plt.ylabel(vert_label)
            if bin_field == 'Density':
                """for density"""
                xlim = plt.xlim()
                plt.xlim(1, xlim[1])
                plt.xscale('log')
                plt.xlabel(r'$\rho$')
                plt.ylabel(vert_label)
            plt.yscale('symlog')
            #plt.legend(loc=0)
            print fname
            plt.savefig(fname)

