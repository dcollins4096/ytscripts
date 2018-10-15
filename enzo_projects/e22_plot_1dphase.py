#!/usr/bin/env python
execfile('go_lite')
import yt
import dsdiff_helpers
#execfile('/Users/dcollins/yt3_scripts/davetools.py')
thelist = ['c21_ct_live',    'c22_ded_live','c20_ppm_live', 'c23_ded_trunk']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 
            'c33_ct_live_field_nodef','c34_ded_live_field_nodef', 'c24_ded_live_def_field_riemann0']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 'c33_ct_live_field_nodef','c34_ded_live_field_nodef',
          'c24_ded_live_def_field_riemann0']
#thelist = ['c34_ded_live_field_nodef', 'c33_ct_live_field_nodef']
thelist = ['c20_ppm_live','c30_ppm_nodef','c37_ppm_def_eta0_Hack2']# 'c36_ppm_def_eta0_defHack']
thelist = ['f01_ppm_nodef', 'f02_ct_nodef', 'f03_ded_nodef']
plot_sequence_name = 'E22_Turb'
base_dir = '/scratch1/dcollins/EnzoProjects/E22_DEF_NonCosmo/'
x_field = 'density'
#y_field = 'magnetic_energy'
#field = 'density'
field = 'kinetic_energy'
field = 'gas_energy'
field = 'ge'
field = 'tvg'
field = 'pe_dave'
field = ("enzo","Density"); ofield = 'enzo_density'
field = ("enzo","TotalEnergy"); ofield = 'enzo_TotalEnergy'
#field = 'TotalEnergy'; ofield = field
#field = 'Temperature'; ofield = 'Temperature'
accumulation = False
def ge_dve(field,data):
    try:
        return data[('enzo','GasEnergy')].v
    except:
        return data[('enzo','TotalEnergy')].v-data['kinetic_energy'].v

def pe_dave(field,data):
    p =data[('enzo','TotalEnergy')].v-0.5*(data[('enzo','x-velocity')].v**2+
                                           data[('enzo','y-velocity')].v**2+
                                           data[('enzo','z-velocity')].v**2)
    p*= data[('enzo','Density')]/(data.ds['Gamma']-1)
    return p
yt.add_field('pe_dave',pe_dave)

yt.add_field('ge',ge_dve)
def temp_vs_ge(field,data):
    return data[('enzo','Temperature')].v/data['ge']

yt.add_field('tvg',temp_vs_ge)

default_frames = [266]
framelist = {'f03_ded_nodef':[125]}
temps = {}
dirnames={}
plt.clf()
for frame in range(3):# range(50): #[0,1,2,5,10,15,20]: #framelist.get(dirname, default_frames):
    plt.clf()
    for dirname in thelist:
        print "frame", frame
        short_name = dirname[0:6]
        setname = dsdiff_helpers.get_ds_name("%s/%s"%(base_dir,dirname),frame)
        print setname
        ds = yt.load(setname)
        #temps[short_name] = yt.create_profile(ds.all_data(),field,fields='cell_mass',weight_field=None)
        prof = yt.create_profile(ds.all_data(),field,fields='cell_mass',weight_field=None,accumulation=accumulation,
                                fractional=True)
        plt.plot(0.5*(prof.x_bins[1:]+prof.x_bins[0:-1]),prof['cell_mass'], label=dirname)
    #plt.xscale('log')
    #plt.yscale('log')
    outname = "%s_%04d_profile_%s"%(plot_sequence_name,frame,ofield)
    if accumulation:
        outname += "cumul"
    plt.legend(loc=3)
    plt.xlabel(ofield)
    plt.ylabel('mass')
    plt.savefig(outname)
    print outname
