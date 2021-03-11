from go import *
reload(taxi)
sims=['p49_1_1','p49_1_half','p49_2_2','p49_3_1','p49_3_half','p49_half_2', 
      'p49_1_2','p49_2_1','p49_2_half','p49_3_2','p49_half_1','p49_half_half']
sims = sims
flt = taxi.fleet(sims)
flt['frames']='all'



#for debugging, everything only uses 0.1% of the domain
#flt['region_type'] = 'rectangle'
#flt['left']  =nar([0.45]*3)
#flt['right'] =nar([0.55]*3)

def add_fluct(obj):
    def bx_fluct(field,data):
        B0 = data.get_field_parameter('B0')
        if  hasattr(B0,'size'):
            B0 = B0[0]
        return data['magnetic_field_x'].v-B0
    obj.add_field('Bx_fluct',function=bx_fluct,sampling_type='cell',validators=[yt.ValidateParameter('B0')])
    def by_fluct(field,data):
        B0 = data.get_field_parameter('B0')
        if  hasattr(B0,'size'):
            B0 = B0[1]
        return data['magnetic_field_y'].v-B0
    obj.add_field('By_fluct',function=by_fluct,sampling_type='cell',validators=[yt.ValidateParameter('B0')])
    def bz_fluct(field,data):
        B0 = data.get_field_parameter('B0')
        if  hasattr(B0,'size'):
            B0 = B0[2]
        return data['magnetic_field_z'].v-B0
    obj.add_field('Bz_fluct',function=bz_fluct,sampling_type='cell',validators=[yt.ValidateParameter('B0')])

if 1:
    for car in flt.taxi_list:
        simname = car.name[4:]
        car.field_parameters={'B0':sim_colors.Means[simname]}

flt['derived_fields']={'fluct':add_fluct}

flt.taxi_list[0].find_extrema(['magnetic_field_z','Bz_fluct'])
print(flt.taxi_list[0].extrema)

if 0:
    allmin=1
    allmax=-1
    for car in flt.taxi_list:
        car.find_extrema(['Bx_fluct','By_fluct','Bz_fluct'])
        print("CAR ====")
        print(car.extrema)
        for fld in ['Bx_fluct','By_fluct','Bz_fluct']:
            allmin = min([allmin, car.extrema[fld][0]])
            allmax = max([allmax, car.extrema[fld][1]])

Bmin = 1e-5
Bmax = 25
bins1 = np.logspace(np.log10(Bmin), np.log10(Bmax),32)
bins2 = -np.logspace(np.log10(Bmin), np.log10(Bmax),32)[::-1]
bins_Bmag = {'magnetic_field_strength':bins1}
#bins_Bx = {'Bx_fluct':np.concatenate([bins2,bins1])}
#bins_By = {'By_fluct':np.concatenate([bins2,bins1])}
#bins_Bz = {'Bz_fluct':np.concatenate([bins2,bins1])}
bins_Bx = {'Bx_fluct':np.linspace(-Bmax,Bmax,256)}
bins_By = {'By_fluct':np.linspace(-Bmax,Bmax,256)}
bins_Bz = {'Bz_fluct':np.linspace(-Bmax,Bmax,256)}

Means = defaultdict(dict)
if 0:
    for car in flt.taxi_list:
        ds = car.load()
        ad = ds.all_data()
        Means[car.name]['Bx'] = ad['magnetic_field_x'].mean()
        Means[car.name]['By'] = ad['magnetic_field_y'].mean()
        Means[car.name]['Bz'] = ad['magnetic_field_z'].mean()


if 1:
    reload(sim_colors)
    for car in flt.taxi_list:
        simname = car.name[4:]
        car.field_parameters={'B0':sim_colors.Means[simname]}

        #car.profile(['magnetic_field_strength','cell_volume'],override_bins=bins_Bmag, field_parameters=field_parameters)

        car.profile(['Bx_fluct','cell_volume'],override_bins=bins_Bx)
        car.profile(['By_fluct','cell_volume'],override_bins=bins_By)
        car.profile(['Bz_fluct','cell_volume'],override_bins=bins_Bz)


