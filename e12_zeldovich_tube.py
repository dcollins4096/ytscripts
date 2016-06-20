if 'ef' not in dir():
    execfile('go')
#ef('tube.py')
import tube
reload(tube)


def _ct_dual_alpha(field,data):
    ke = 0.5*(data[('enzo','x-velocity')]**2)
    try: 
        ke += 0.5*data[('enzo','y-velocity')]**2 
    except:
        pass
    try:
        ke += 0.5*data[('enzo','z-velocity')]**2 
    except:
        pass
    try:
        be = 0.5*(data[('enzo','Bx')]**2+data[('enzo','By')]**2 + data[('enzo','Bz')]**2)/data[('enzo','Density')]
        be=be.v
    except:
        be = 0.0
    all_ke = ke.v + be
    thermal  = (data[('enzo','TotalEnergy')].v - all_ke )
    alpha = thermal/all_ke/0.008
    return alpha
yt.add_field('ct_dual_alpha',_ct_dual_alpha, take_log=False)
def _enzo_pressure(field,data):
    ke = 0.5*(data[('enzo','x-velocity')]**2)
    try: 
        ke += 0.5*data[('enzo','y-velocity')]**2 
    except:
        pass
    try:
        ke += 0.5*data[('enzo','z-velocity')]**2 
    except:
        pass
    try:
        be = 0.5*(data[('enzo','Bx')]**2+data[('enzo','By')]**2 + data[('enzo','Bz')]**2)/data[('enzo','Density')]
        be=be.v
    except:
        be = 0.0
    p = (data[('enzo','TotalEnergy')].v -ke.v - be)*(data.ds['Gamma']-1)*data[('enzo','Density')].v
    return p
yt.add_field('enzo_pressure',_enzo_pressure,take_log=False)
def _plasma_beta(field,data):
    try:
        return data['pressure'].v/(0.5*data['magnetic_energy'].v)
    except:
        return data['pressure']*0
#yt.add_field('plasma_beta',function=_plasma_beta)
def _BoverRho(field,data):
    b_norm = np.sqrt((data[('enzo','Bx')]**2+data[('enzo','By')]**2 + data[('enzo','Bz')]**2)).v
    return b_norm/data[('enzo','Density')].v
#yt.add_field('BoverRho',function=_BoverRho)
def _PoverRho(field,data):
    return data['enzo_pressure']/data[('enzo','Density')].v
yt.add_field('PoverRho',_PoverRho)
def _sound_squared(field,data):
    return data.ds['Gamma']*data['PoverRho']
yt.add_field('cs2',_sound_squared)
def _UseDualEnergy(field,data):
    out = -1*np.ones_like(data['cs2']) - 0.01
    eta1 = 1e-3
    ke = 0.5*(data[('enzo','x-velocity')]**2)
    try: 
        ke += 0.5*data[('enzo','y-velocity')]**2 
    except:
        pass
    try:
        ke += 0.5*data[('enzo','z-velocity')]**2 
    except:
        pass
    try:
        be = 0.5*(data[('enzo','Bx')]**2+data[('enzo','By')]**2 + data[('enzo','Bz')]**2)/data[('enzo','Density')]
        be=be.v
    except:
        be=0.0
    UseDEF = np.logical_and(data['cs2'] > eta1*2*ke.v, data['cs2'] > eta1*2.*be) #2 is for v^2 vs. ke.
    #There is another conditional that comes from the solver, and not accessible from here:
    # int_from_gas_solver < 0.5*int_from_ke
    out[UseDEF] = -0.1
    return out
yt.add_field('UseDualEnergy',_UseDualEnergy,take_log=False)


z01 = 'z01_ppm_fiducial'
z02 = 'z02_ded_fiducial'
z03 = 'z03_ct_fiducial'
z01b = 'z01b_ppm_nodef'
z02b = 'z02b_ded_nodef'
z03b = 'z03b_ct_nodef'
z04 = 'z04_ded_smallb'
z04b = 'z04_ded_smallb_backup_oldexe'
z05 = 'z05_ct_smallb'
z05b = 'z05_ct_smallb_backup'
z06 = 'z06_ded_largeb'
z07 = 'z07_ct_largeb'
z08 = 'z08_ded_bw'
z09 = 'z09_ct_bw'
z10b = 'z10b_ded_horse'
z11 = ''
z12 = 'z12_ded_collapse'
z13 = 'z13_ct_collapse'
z14 = 'z14_ded_coll_def'
z15 = 'z15_ct_coll_def'
z16 = 'z16_ded_coll_nodef_smallp'
z17 = 'z17_ct_coll_nodef_smallp'
z18 = 'z18_ded_coll_def_smallp'
z19 = 'z19_ct_col_def_smallp'
z20 = ''
z21 = 'z21_ct_shock'
z30 = 'z30_ded_largeb_old_dbda'
z31 = 'z31_ded_largeb_expansion1'
z32 = 'z32_ct_largeb_old_deda'
z33 = 'z33_ct_largeb_new_deda'
z40 = 'z40_ded_largefield_novelocity'
z41 = 'z41_ct_largefield_novelocity'
z42 = 'z42_ded_largefield_hightemp'
z43 = 'z43_ct_largefield_hightemp'

#z06 = 'z06_ded_hack' #the way I think it should be
#z07 = 'z07_ct_hack'  #ct
#z08 = 'z08_ded_no_hack' #the way it is
#z09 = 'z09_ded_energy_only'
#z10 = 'z10_ded_goofoff'
#z11 = 'z11_ct_goofoff'
#frames = range(22)
#frames = [21]
#frames = [0,11]
frames = [100]
#frames = range(0,50,5)

fields = ['Density','x-velocity', 'TotalEnergy', 'Temperature','enzo_pressure'] #,'By']
#my_sims = [z01, z04,z05]; note = 'small_b' #_with_b'
#my_sims = [z04,z05]; note = 'small_b_with_b'; fields += ['By']
#my_sims = [z02,z03,z04,z05]; note = 'z02-z05-withb'; fields += ['By']
#my_sims = [z04,z04b]; note = 'test'; fields += ['By']
#my_sims = [z05,z05b]; note = 'tes5t'; fields += ['By']
#my_sims = [z01,z02,z03]
#my_sims = [z04b,z06]; note = 'version_test'; fields += ['By']
#my_sims = [z06,z07]; note = 'large_b_old_expansion_temp'; fields += ['By', 'GasEnergy']
#my_sims = [z08,z09]; note = 'brio'; fields += ['By', 'PoverRho']
#my_sims = [z08,z09]; note = 'brio_dual'; fields += ['By', 'GasEnergy']
#my_sims = [z10b]; note = 'collapse' 
#my_sims = [z12,z13]; note = 'collapse'; fields += ['By'] #,'GasEnergy']
#my_sims = [z14,z15]; note = 'collapse_def'; fields += ['By','GasEnergy']
#my_sims = [z16,z17]; note = 'coll_nodef_smallp'; fields += ['By']; frames =  range(100)
#my_sims = [z18,z19]; note = 'coll_def_smallp'; fields += ['By', 'GasEnergy']; frames = [100] #  range(100)
#my_sims = [z21]; note = 'shock'; fields += ['By','ct_dual_alpha']; frames = [1] #  range(100)

#Use this set to convince self that the new dbda is wrong (the old version is correct.)
#my_sims = [z07,z30,z31,z32,z33]; note = 'dbda_test' ; fields += ['By']; frames=[21]

#Testing for the PR:  Why does pressure become negative in CT?
#my_sims = [z30,z32]; note = 'DEF_test' ; fields += ['By', 'UseDualEnergy']; frames=[0] #range(22)
#my_sims = [z30,z32]; note = 'DEF_test_ponly' ; fields = ['enzo_pressure']; frames=range(22)
#These look good, pressure non-negative
#my_sims = [z01,z02,z03]; note = 'DEF_test_ponly_nofield' ; frames=[21]; #fields = ['enzo_pressure']; frames=[21]
#my_sims = [z01,z04,z05]; note = 'DEF_test_ponly_smallfield' ; frames=[21]; #fields = ['enzo_pressure']; frames=[21]
#my_sims = [z40,z41]; note = 'noV' ; frames=[0,21]; #fields = ['enzo_pressure']; frames=[21]
#high temp sims look good.
my_sims = [z42,z43]; note = 'arg' ;  fields += ['By'] ; 
frames = [0,1]

basedir = '/Users/dcollins/scratch/EnzoProjects/E12/E12c_Zeldovich'
sim_list = ["%s/%s"%(basedir,s1) for s1 in my_sims]
labels = [s1[0:7] for s1 in my_sims]
delta = False
#delta = True; note += 'delta'

OutputTemplate = "/DD%04d/data%04d"
#OutputTemplate = "/RD%04d/RedshiftOutput%04d"

axis=1
counter=0
#fields=['density','TotalEnergy','x-velocity','By','enzo_pressure']
#renorm={'density':1,'pressure':0.6,'TotalEnergy':0.9}

for frame in frames:

    ds_list = [ yt.load(sim_1+OutputTemplate%(frame,frame)) for sim_1 in sim_list ]

    y=tube.tube(ds_list,  delta=delta, fields=fields, renorm=False,
                filename = 'e12_tube_%s_cons_%04d.pdf'%(note,frame), legend=True, labels=labels, xlim=[0.42,0.58])
if 0:
    fields_derived=['thermal_energy', 'pressure','BoverRho']
    y=tube.tube(ds_list,  delta=delta, fields=fields_derived, renorm=False,
                filename = 'e12_tube_%s_derived_%04d.pdf'%(note,frame), legend=True, labels=labels)
