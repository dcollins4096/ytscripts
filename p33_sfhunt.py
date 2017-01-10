if 'ef' not in dir():
    execfile('go')
reload(taxi)

def _ref_mass(field,data):
    """
  float ModifiedMinimumMassForRefinement =
    MinimumMassForRefinement[method]*POW(RefineBy,
            level*MinimumMassForRefinementLevelExponent[method]);
	FlaggingField[i] = ((ffield[i] > ModifiedMinimumMassForRefinement)
	(then some other stuff that I don't care about right now.)
	"""
    method = np.where(data.ds['CellFlaggingMethod'] == 2)[0][0]
    exp = data.ds['MinimumMassForRefinementLevelExponent'][method]
    ModifiedMinimumMassForRefinement=data.ds['MinimumMassForRefinement'][method]*2**(data.Level*exp)
    output = np.zeros(data['cell_mass'].shape)
    flag = data['cell_mass'].in_units('code_mass') > ModifiedMinimumMassForRefinement
    output[flag] = 1
    return output
yt.add_field('ref_mass',function=_ref_mass,take_log=False,validators=[yt.ValidateGridType()])

def _ref_metal(field,data):
    output = np.zeros(data['Metal_Density'].shape)
    flag = data['Metal_Density']/data['density'].in_units('code_density')/0.022 > data.ds['MetallicityRefinementMinMetallicity']
    output[flag]=1
    return output
yt.add_field('ref_metal',function=_ref_metal,take_log=False,validators=[yt.ValidateGridType()])

    

def _sf_finest_level(field,data):
    output = np.zeros(data['density'].shape)
    if hasattr(data,'child_mask'):
        output[data.child_mask] = 1
    return output
yt.add_field('sf_finest_level',function=_sf_finest_level,take_log=False, validators=[yt.ValidateGridType()])

def _sf_divv(field,data):
    output = np.zeros(data['velocity_divergence'].shape)
    output[ data['velocity_divergence'] < 0 ] = 1
    return output
yt.add_field('sf_divv',function=_sf_divv,take_log=False)

def _sf_overdensity(field,data):
    output = np.zeros(data['density'].shape)
    odthresh = data.ds['StarMakerOverDensityThreshold']
    output[ data['density'] > odthresh ] = 1
    return output
yt.add_field('sf_overdensity',function=_sf_overdensity,take_log=False)

def _sf_timescale(field,data):
    G = data.ds.quan(6.67428e-8, '1/(code_density*s**2)')
    d1 = data.ds['DensityUnits']
    t1 = data.ds['TimeUnits']
    d = data['density'].in_units('code_density')*d1 #might be round-about, but it ensures similarity
    tdyn = np.sqrt(3*np.pi/(32*G*d))/t1
    output = np.ones(data['density'].shape)
    nostar = np.logical_and( tdyn < data['enzo','Cooling_Time'], data[('enzo','Temperature')] > 1e4)
    output[nostar]=0
    return output
yt.add_field('sf_timescale',function=_sf_timescale,take_log=False)

def _sf_jeans(field,data):
    """
               bmass = d(i,j,k)*dble(d1)*dble(x1*dx)**3 / msolar
               isosndsp2 = sndspdC * temp(i,j,k)
               jeanmass = pi/(6._RKIND*sqrt(d(i,j,k)*dble(d1))) *
     &                    dble(pi * isosndsp2 / G)**1.5_RKIND / msolar

               if (bmass .lt. jeanmass) goto 10z
	"""
    msolar = data.ds.quan(1.9891e33,'g')
    G = data.ds.quan(6.67428e-8, 'cm**3/(g*s**2)')
    d1 = data.ds['DensityUnits']
    t1 = data.ds['TimeUnits']
    L1 = data.ds['LengthUnits']
    d = data['density'].in_units('g/cm**3') #data['density'].in_units('code_density')*d1 #might be round-about, but it ensures similarity
#bmas=d*(data['dx'][0,0,0].in_units('code_length')*L1)**3/msolar
    bmas=d*(data['dx'][0,0,0].in_units('cm'))**3/msolar
    T = data[('enzo','Temperature')]
    if hasattr(T,'units'):
        k_over_m_units = 'cm**2/s**2/K'
    else:
        k_over_m_units = 'cm**2/s**2'
    sndspdC=data.ds.quan(1.3095e8, k_over_m_units) #kB/(0.6 proton mass)
    isosndsp2 = sndspdC*T
    jeansmass = np.pi/(6*np.sqrt(d))*(np.pi*isosndsp2/G)**1.5/msolar
    output = np.zeros(bmas.shape)
    output[bmas>jeansmass] = 1
    return output
yt.add_field('sf_jeans',function=_sf_jeans,take_log=False,validators=[yt.ValidateGridType()])

if 1:
    aj15=taxi.taxi('aj15_sphere')
    aj16=taxi.taxi('aj16_sphere')
    aj17=taxi.taxi('aj17_sphere')
    fleet  = taxi.fleet([aj15,aj16,aj17])
    fleet['weight_field']='dx'

all_sf_fields =['sf_finest_level','sf_divv','sf_overdensity','sf_timescale','sf_jeans']
all_ref_fields = ['ref_mass','ref_metal']
fleet['fields']=all_ref_fields
fleet['frames']=[0,1]
fleet.plot()
#print aj15.ds.index.grids[0]
#reg = aj15.get_region()
#stat(reg['sf_divv'])
