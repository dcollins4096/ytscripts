#i22=taxi.taxi('i22')
i22.fields=['Bstretching_f2', 'Bstretching_p2','Bstretching_m2']
i22.fields+=['Badvection_f2', 'Badvection_p2','Badvection_m2']
i22.fields+=['Bcompression_f2', 'Bcompression_p2','Bcompression_m2']
for field in i22.fields:
    print i22.ds.index.grids[-1][field][0]
#i22.fields=[ 'Bstretching_p2','Bstretching_m2']
i22.frames=[70]
i22.plot()

if 'fleet' not in dir() and False:
    fleet = taxi.fleet(['i10','i02'])
    avg_b2={}
    max_b2={}

    t1 = {}
    npart = {}
    ncycle={}
    mass={}
    frames={}

if 0:
    for car in fleet.taxi_list:
        t1[car.name] = t1.get(car.name, [])
        npart[car.name] = npart.get(car.name,[])
        ncycle[car.name]= ncycle.get(car.name,[])
        mass[car.name]= mass.get(car.name,[])
        frames[car.name] = frames.get(car.name,[])
        avg_b2[car.name] = avg_b2.get(car.name,[])
        max_b2[car.name] = max_b2.get(car.name,[])
        for n in car.frames:
            if n in ncycle[car.name]:
                print "SKIPPING", n, ncycle, car.name
                continue
            car.fill(n)
            ncycle[car.name].append(car.ds['InitialCycleNumber'])
            t1[car.name].append(car.ds['InitialTime'])
            ad=car.ds.all_data()
            avg_b2[car.name].append(ad.quantities['WeightedAverageQuantity']('magnetic_field_strength','cell_volume'))
            max_b2[car.name].append(ad['magnetic_field_strength'].max())


if 1:
    carnames = "all_%s"%fleet.allnames()
    plt.clf()
    colordict = {'i02':'b','i10':'g'}
    for car in fleet.taxi_list:
        plt.plot(t1[car.name],avg_b2[car.name],label=r"$\bar{b}$ %s"%car.outname,marker='x',c=colordict.get(car.name,'k'))
        plt.plot(t1[car.name],max_b2[car.name],label=r"$b_{max}$ %s"%car.outname,marker='x',linestyle='--',c=colordict.get(car.name,'k'))
    plt.legend(loc=0)
    plt.xlabel('time(code)'); plt.ylabel('B')
    plt.yscale('log')
    outname = 'p33_%s_field_time.pdf'%carnames
    plt.savefig(outname)
    print outname



