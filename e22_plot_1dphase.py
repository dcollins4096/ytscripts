#!/usr/bin/env python
execfile('/Users/dcollins/yt3_scripts/go_lite')
import yt
#execfile('/Users/dcollins/yt3_scripts/davetools.py')
basedir ='/Users/dcollins/scratch/EnzoProjects/E12/E12b_Cosmology/' 
thelist = ['c21_ct_live',    'c22_ded_live','c20_ppm_live', 'c23_ded_trunk']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 
            'c33_ct_live_field_nodef','c34_ded_live_field_nodef', 'c24_ded_live_def_field_riemann0']
thelist += ['c30_ppm_nodef', 'c31_ct_live_nofield', 'c32_ded_live_nofield', 'c33_ct_live_field_nodef','c34_ded_live_field_nodef',
          'c24_ded_live_def_field_riemann0']
#thelist = ['c34_ded_live_field_nodef', 'c33_ct_live_field_nodef']
#thelist = ['c20_ppm_live','c30_ppm_nodef','c37_ppm_def_eta0_Hack2']# 'c36_ppm_def_eta0_defHack']
#thelist = ['c01b_ppm_hll','c02b_ct_hll','c03b_ded_hll']
#thelist = ['c60_ppm_allthings',       'c61_ct_allthings',        'c62_ded_allthings']#,'c02_ct', 'c01_ppm']
#thelist = ['c64_ded_bsmall','c63_ct_bsmall','c65_ct2_bsmall']
#thelist = ['c01_ppm'] #,'c60_ppm_allthings','c61_ct_allthings', 'c02_ct']
#thelist += ['c66_ded_bsmall_Rad12','c67_ct2_bsmall_rad12','c68_ppm_rad12'] #RadiationFieldType = 12
#thelist += ['c61_ct_allthings','c62_ded_allthings']
#thelist = ['c01_ppm'] #,'c60_ppm_allthings','c61_ct_allthings', 'c02_ct']
thelist = ['c66_ded_bsmall_Rad12','c67_ct2_bsmall_rad12','c68_ppm_rad12'] #RadiationFieldType = 12
thelist = ['c69_hydro3_rad12','c66_ded_bsmall_Rad12','c67_ct2_bsmall_rad12','c68_ppm_rad12'] #RadiationFieldType = 12
plot_sequence_name = 'Enzo15_1d'
x_field = 'density'
#y_field = 'magnetic_energy'
#field = 'density'
field = 'kinetic_energy'
field = 'gas_energy'
field = 'ge'
field = 'tvg'
field = 'pe_dave'
ofield = None
#field = ("enzo","Density"); ofield = 'enzo_density'
#field = ("enzo","TotalEnergy"); ofield = 'enzo_TotalEnergy'
field = 'TotalEnergy'
field = 'Temperature'
#field = 'magnetic_energy'
#field = 'Electron_Density'
field = 'HI_Density'
if 'field_loop' in dir():
    field = field_loop
    #print "CLOWN ok"
    #pdb.set_trace()

#['HII_Density','HI_Density','HeIII_Density','HeII_Density','HeI_Density','Metal_Density']

outname_trail = ''
#field = 'density'

if ofield is None:
    ofield = field
accumulation = False
print field

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

temps = {}
dirnames={}
plt.clf()

print "CLOWN",field
offset = 1.0
for frame in [0]:
    print "frame", frame
    for dirname in thelist:
        short_name = dirname[0:6]
        outname_trail += "_%s"%short_name
        print "prof", short_name
        #ds = yt.load("%s/%s/DD%04d/DD%04d"%(basedir,dirname,frame,frame))
        ds = yt.load("%s/%s/RD%04d/RD%04d"%(basedir,dirname,frame,frame))
        #temps[short_name] = yt.create_profile(ds.all_data(),field,fields='cell_mass',weight_field=None)
        prof = yt.create_profile(ds.all_data(),field,fields='cell_mass',weight_field=None,accumulation=accumulation,
                                fractional=True)
        plt.plot(0.5*(prof.x_bins[1:]+prof.x_bins[0:-1]),offset*prof['cell_mass'], label=dirname,linewidth=1)
        #offset += 0.1
    plt.xscale('log')
    plt.yscale('log')
    outname = "%s_%04d_profile_%s%s.pdf"%(plot_sequence_name,frame,ofield,outname_trail)
    if accumulation:
        outname += "cumul"
    plt.legend(loc=4)
    plt.xlabel(ofield)
    plt.ylabel('mass')
    plt.savefig(outname)
    print outname
