
#!/usr/bin/env python
execfile('/Users/dcollins/yt3_scripts/go_lite')
import yt
#execfile('/Users/dcollins/yt3_scripts/davetools.py')
def stat(array,strin='', format='%0.16e'):
    template = '['+format+','+format+'] %s %s'
    print template%(array.min(),array.max(),array.shape,strin)

thelist = ['c21_ct_live',    'c22_ded_live','c20_ppm_live', 'c23_ded_trunk']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 
            'c33_ct_live_field_nodef','c34_ded_live_field_nodef', 'c24_ded_live_def_field_riemann0']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 'c33_ct_live_field_nodef','c34_ded_live_field_nodef',
          'c24_ded_live_def_field_riemann0']
#thelist = ['c34_ded_live_field_nodef', 'c33_ct_live_field_nodef']
thelist = ['c20_ppm_live','c30_ppm_nodef']
thelist = ['c20_ppm_live','c30_ppm_nodef','c37_ppm_def_eta0_Hack2']# 'c36_ppm_def_eta0_defHack']
#thelist += ['c38_ppm_def_eta0_hack3']
#thelist += ['c39_ppm_def_eta0_DEFhack4']
#thelist = ['x02_nodef_unigrid', 'x01_def_unigrid','x04_ELhaxk','x05_hax']
#thelist = ['x02_nodef_unigrid','x04_ELhaxk','x06_hax_more'] #'x05_hax']
#thelist = ['x01_def_unigrid', 'x02_nodef_unigrid']
thelist = ['x09_def_etaDefault_noPmin', 'x02_nodef_unigrid']
plot_sequence_name = 'PPM_etaDefault_b_nopmin'
x_field = 'density'
#y_field = 'magnetic_energy'
#field = 'density'
field = 'TotalEnergy'
field = 'kinetic_energy'
field = 'gas_energy'
field = 'ge'
field = 'Temperature'
lots = True
accumulation = False


for frame in range(5)+[5,10, 100, 200, 278]: #range(10):
    if lots:
        tmp_fig = plt.figure()
        tmp_ax = tmp_fig.add_subplot(111)
        ke_fig = plt.figure()
        ke_ax = ke_fig.add_subplot(111)
        te_fig = plt.figure()
        te_ax = te_fig.add_subplot(111)
        ge_fig2 = plt.figure()
        ge_ax2 = ge_fig2.add_subplot(111)
    for dirname in thelist:
        short_name = dirname[0:6]
        print "prof", short_name
        #ds = yt.load()
        fname = "%s/DD%04d/DD%04d.cpu0000"%(dirname,frame,frame)
        fptr = h5py.File(fname)
        try:
#        if True:
            grid = fptr['Grid%08d'%(1)]
            temperature = grid['Temperature'][:]
            ge_disk = None
            if lots:
                density = grid['Density'][:]
                vx = grid['x-velocity'][:]
                vy = grid['y-velocity'][:]
                vz = grid['z-velocity'][:]
                te = grid['TotalEnergy'][:]
                if 'GasEnergy' in grid.keys(): 
                    print "I have gas"
                    ge_disk = grid['GasEnergy'][:]
                    stat(ge_disk,"Gas Energy, Disk")
                else:
                    print "Gas-X.  Its for your butt."
                ke = 0.5*(vx*vx+vy*vy+vz*vz)
                ge_find = te - ke
                stat(ge_find,"Gas energy, Derived.")

        except:
            raise
        finally:
            pass
            fptr.close()
        bins = 100
        bins = np.arange(0, 7,7./100)
        tmp_ax.hist(np.log10(temperature.flatten()),histtype='step',label=dirname, bins=bins)
        tmp_ax.set_yscale('log'); tmp_ax.set_xlabel('log temp disk'); tmp_ax.legend(loc=2)
        if lots:
            ke_ax.hist(np.log10(ke.flatten()), histtype='step',label=dirname, bins=100)
            ke_ax.set_xlabel('log10 ke')
            ke_ax.set_yscale('log'); ke_ax.legend(loc=2)
            te_ax.hist(np.log10(te.flatten()), histtype='step',label=dirname, bins=100)
            te_ax.set_xlabel('log10 te')
            te_ax.set_yscale('log'); te_ax.legend(loc=2)
            ge_ax2.hist((ge_find.flatten()), histtype='step',label=dirname, bins=100)
            ge_ax2.set_xlabel('log10 ge derived')
            if ge_disk is not None:
                ge_ax2.hist(ge_disk.flatten(), histtype='step',label=dirname+' dsk', bins=100)
            ge_ax2.set_yscale('log'); ge_ax2.legend(loc=2)
    outname = "%s_%04d_disk"%(plot_sequence_name,frame)

    if lots:
        ge_fig2.savefig(outname + "_GE")
        ke_fig.savefig(outname + "_KE")
        te_fig.savefig(outname + "_TE")

    if 1:
        tmp_fig.savefig(outname+"_Temp")
        print outname



