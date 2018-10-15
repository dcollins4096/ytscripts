execfile('go')
execfile('go')  
mhd=True            
ef('../yt3_scripts/p12_fields_rectified.py')
ef('../yt2.0/p12_Laplacians.py')
for simset in [0]: #[0,1]:
    
    if simset==1:       
        frame = 60; basedir = '/scratch1/dcollins/Paper08/B02/512/'
        peaks = fPickle.load('../yt3_scripts/p12_b02_512_n0060_peaks2.pickle'); sim='b025v'
    if simset==0:   
        frame = 60; basedir = '/scratch1/dcollins/Paper08/B20/512/'
        peaks = fPickle.load('../yt3_scripts/p12_b20_512_n0060_peaks2.pickle'); sim='b205v'
                    
    fname = '%s/RS%04d/restart%04d'%(basedir, frame,frame)
    t0 = time.time()

    ds = yt.load(fname)
#    if 'mhd' not in dir():
#        ef('p13_fields.py')
    fields = ['Lambda_therm_pos', 'Lambda_therm_neg', 'Lambda_turb_pos','Lambda_turb_neg',
              'Lambda_magn_pos', 'Lambda_magn_neg']
    Colors = {'Lambda_therm_pos':'b', 'Lambda_therm_neg':'b', 'Lambda_turb_pos':'r','Lambda_turb_neg':'r',
              'Lambda_magn_pos':'g', 'Lambda_magn_neg':'g'}
    Lines  = {'Lambda_therm_pos':'-', 'Lambda_therm_neg':'--', 'Lambda_turb_pos':'-','Lambda_turb_neg':'--',
              'Lambda_magn_pos':'-', 'Lambda_magn_neg':'--'}
    AllPlots = {}
    for field in fields:
        AllPlots[field]=[]

    all_fields = ["CellVolumeNorm", "Delta", "Gravity", "KineticEnergyDensity", "Pressure", 
           "GradProductDensityPressure", "LaplacePressure", "Lambda_therm", "Lambda_therm_pos", 
           "Lambda_therm_neg", "EnstrophyDensity", "RateOfStrainSqr", "Denstrophy", "DiffOmegaSqrStrainSqr", 
           "Lambda_turb", "Lambda_turb_pos", "Lambda_turb_neg", "MagneticPressure", "GradProductDensityMagneticPressure", 
           "LaplaceMagneticPressure", "MagneticFieldCrossContractions", "Lambda_magn", "Lambda_magn_pos", "Lambda_magn_neg"]

    t0 = time.time()


#   fields = ['Lambda_therm','Lambda_turb','Lambda_magn','density','velocity_divergence']
#   fields = ['Lambda_therm_per_delta','Lambda_turb_per_delta','Lambda_magn_per_delta','Lambda_per_delta']
#   fields=['Gravity']
#   fields=['Lambda_therm']
#   fields=['GradProductDensityPressure','LaplacePressure']
#   fields=['partial_lambda_grav']
#   fields=['Gravity']
#   fields=['Lambda_magn','LaplaceMagneticPressure',
#           'GradProductDensityMagneticPressure','MagneticFieldCrossContractions']
#   fields =['LaplaceMagneticPressure_x','LaplaceMagneticPressure_y','LaplaceMagneticPressure_z']
    fields=['Lambda_therm_full_per_G_pos','Lambda_therm_full_per_G_neg' ,
                              'Lambda_magn_per_G_pos','Lambda_magn_per_G_neg',
                              'Lambda_turb_per_G_pos','Lambda_turb_per_G_neg','density','velocity_divergence']
    #fields=['LaplacianMagneticPressure_27pt']
    fields=['Lambda_magn']
    Fcol = {'Lambda_therm':'b','Lambda_turb':'r','Lambda_magn':'g','density':'k','velocity_divergence':'c'}
    Fcol.update({'Lambda_therm_pos_per_delta':'b', 'Lambda_turb_pos_per_delta':'r','Lambda_turb_neg_per_delta':'r',
                 'Lambda_therm_neg_per_delta':'b','Lambda_magn_pos_per_delta':'g','Lambda_magn_neg_per_delta':'g',
                 'Lambda_therm_per_delta':'b','Lambda_magn_per_delta':'g','Lambda_turb_per_delta':'r'})
    Line = {'Lambda_turb_neg_per_delta':'--', 'Lambda_therm_neg_per_delta':'--','Lambda_magn_neg_per_delta':'--',
            'Lambda_turb_pos_per_delta':'-.', 'Lambda_therm_pos_per_delta':'-.','Lambda_magn_pos_per_delta':'-.'} 


    symlog_thresh = {'Lambda_therm':1,'Lambda_turb':1,'Lambda_magn':1,'velocity_divergence':1,
                     'Lambda_magn':1,'LaplaceMagneticPressure':1, 'GradProductDensityMagneticPressure':1,
                      'MagneticFieldCrossContractions':1}
    labels = dict(zip(['Lambda_therm','Lambda_turb','Lambda_magn','density','velocity_divergence'],
                      ['LH','LT','LM','rho','d']))
    Fcol['Lambda_therm_per_delta']='b'
    labels['Lambda_therm_per_delta'] = 'LH/d'
    labels.update({'Lambda_therm_pos_per_delta':'LHp/d', 'Lambda_turb_pos_per_delta':'LTp/d','Lambda_turb_neg_per_delta':'LTn/d',
                 'Lambda_therm_neg_per_delta':'LHn/d','Lambda_magn_pos_per_delta':'LMp/d','Lambda_magn_neg_per_delta':'LMp/d',
                 'Lambda_therm_per_delta':'LHt/d','Lambda_magn_per_delta':'LMt/d','Lambda_turb_per_delta':'LTt/d'})
    Rlist = [0.5e-3, 1e-3, 1e-2, 5e-2]
    Rcolo = [   'r',  'k',  'r',  'k']
    ax=0
    symlog_thresh['LaplaceMagneticPressure_27pt_bh'] = 1
    for npeak, peak in enumerate(peaks):
        if npeak != 0:
            continue
        width = 0.05
        data = ds.h.sphere(peak,(width,'code_length'))

        for field in fields:
            proj = ds.proj(field, ax, data_source = data, center = peak)
            pw = proj.to_pw(center = peak, width = 2*width)
            for ncol,R in enumerate(Rlist):
                pw.annotate_sphere(peak, (R,'code_length'), circle_args={'color':Rcolo[ncol]})
            if symlog_thresh.has_key(field):
                pw.set_log(field,True,linthresh=symlog_thresh[field])
            if field not in ['density']:
                #pw.set_log(field,False)
                #pw.set_zlim(field,-5,5)
                pw.save('%s_n%04d_peak%04d'%(sim,frame,npeak))
            else:

                pw.annotate_magnetic_field()
                pw.save('b025_n%04d_peak%04d_MAGNETIC'%(frame,npeak))
                pw._callbacks=[]
                for ncol,R in enumerate(Rlist):
                    pw.annotate_sphere(peak, (R,'code_length'), circle_args={'color':Rcolo[ncol]})
                pw.annotate_velocity()
                pw.save('%s_n%04d_peak%04d_VELOCITY'%(sim,frame,npeak))

                    
