if 0:
    phase_args={'y_bounds':nar([1e-17,4e-12]),'x_bounds':nar([1e-2,1e4])}
    phase_args['x_bins']=128
    phase_args['y_bins']=128
    reload(davecallback)
    cb = [davecallback.Equlib()]
    cb += [davecallback.ConstTempLines()]
    oober.phase(['davesNumberDensity','GasPressure','CellMass'],phase_args=phase_args,callbacks=cb)
    #oober.pc.set_zlim(1e39,1e40)
    #print oober.pc.save('arse')
    #dbounds = oober.region.quantities['Extrema']('davesNumberDensity')[0]
    #pbounds = oober.region.quantities['Extrema']('GasPressure')[0]
    #if dbounds[0] < phase_args['x_bounds'][0] or dbounds[1] > phase_args['x_bounds'][1]:
        #print "d bounds", dbounds
    #if pbounds[0] < phase_args['y_bounds'][0] or pbounds[1] > phase_args['y_bounds'][1]:
        #print "p bounds", pbounds
    #
    #cool4.pc.plots[0]._callbacks=[cb]
    #print cool4.pc.save('arse')

if 1:
    phase_args={'y_bounds':nar([10**(-1),10**(5.5)]),'x_bounds':nar([1e-2,1e4])}
    phase_args['x_bins']=128
    phase_args['y_bins']=128
    reload(davecallback)
    cb = [davecallback.Equlib()]
    cb += [davecallback.ConstTempLines(temp_list=[10,18,100,1000,1e4,5250])]

if 1:
    #cb += [davecallback.title(title='Central octant')]
    oober.plot()
    oober.phase(['davesNumberDensity','PoverK','CellMass'],phase_args=phase_args,callbacks=cb)
    #oober.phase(['davesNumberDensity','PoverK','CellVolume'],phase_args=phase_args,callbacks=cb)

if 0:
    save = oober.outname
    oober.outname += 'VolumeWeight'
    oober.phase(['davesNumberDensity','PoverK','LorentzForceNorm'],
                weight='CellVolume',phase_args=phase_args,callbacks=cb)
    oober.outname = save
    save = oober.outname
    oober.outname += 'Mass'
    oober.phase(['davesNumberDensity','PoverK','LorentzForceNorm'],
                weight='CellMassMsun',phase_args=phase_args,callbacks=cb)
    oober.outname = save

if 0:
    oober.phase(['daveTemperature','LorentzAccel','CellMass'])
if 0:
    x1 = plot_phase('ProfileFiles/MP2_t9_0001_daveTemperature_LorentzAccel_CellMass_0.pickle')
    x12= plot_phase('ProfileFiles/MP2_t9_0012_daveTemperature_LorentzAccel_CellMass_0.pickle')
    for x in (x1,x12):
        x[1].set_xlim(10,1e5)
        x[1].set_ylim(1e-11,1e-4)
    x1[1].save_image('MP2_t9_0001_Profile2D_0_daveTemperature_LorentzAccel_CellMass')
    x12[1].save_image('MP2_t9_0012_Profile2D_0_daveTemperature_LorentzAccel_CellMass')
        

if 0:
    #phase_args={'y_bounds':nar([10**(2.5),10**(5.5)]),'x_bounds':nar([1e-2,1e4])}
    phase_args={'y_log':False}
    phase_args['x_bins']=128
    phase_args['y_bins']=128
    reload(davecallback)
    cb = [davecallback.Equlib()]
    cb += [davecallback.ConstTempLines()]
    oober.phase(['davesNumberDensity','PoverK','CellMass'],phase_args=phase_args,callbacks=cb)
