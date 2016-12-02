
aq32 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq32_m2.9_drive1_32',name='aq32',frames=range(0,42,2),fields=['density'])
aq31 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq31_m2.9_drive0.5_32',name='aq31',frames=range(0,42,2),fields=['density'])
aq16 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq16_ppm_m9_drive0_noamr_128',name='aq15',frames=range(17),fields=['density'])
aq33 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq33_m9_drive0.5_32',name='aq33',frames=range(0,110,10))
aq34 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq34_m9_drive0_p60_32', name='aq34', frames=range(0,110,10))
aq38 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq38_m2.9_drive0_32_p59_ppm',name='aq38',frames=range(0,42,2),fields=['density'])
aq39 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq39_m2.9_drive0_32_p59_ppm' ,name='aq39',frames=range(0,42,2),fields=['density'])
aq40 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq40_m2.9_drive0_32_p59_hydro6_newAddVel' ,name='aq40',frames=range(0,42,2),fields=['density'])
aq41 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq41_m9_drive0.5_32_p59_ppm' ,name='aq41',frames=range(0,42,2),fields=['density'])


fleet = [aq41]
for car in fleet:
    car.frames = range(41)
    car.plot()
    print car.name
    vx=[]
    vy=[]
    vz=[]
    mach=[]
    px=[]
    py=[]
    pz=[]
    ex=[]
    ey=[]
    ez=[]
    t = []
    tdyn = 0.5/9
    
    for frame in car.frames:
        car.fill(frame)
        ds=car.ds
        t.append(car.ds.current_time/(tdyn) )
        reg = car.get_region(frame)
        vx.append(reg.quantities['WeightedAverageQuantity']('velocity_x','cell_volume'))
        vy.append(reg.quantities['WeightedAverageQuantity']('velocity_y','cell_volume'))
        vz.append(reg.quantities['WeightedAverageQuantity']('velocity_z','cell_volume'))
        mach.append( np.sqrt(reg.quantities['WeightedAverageQuantity']('mean_square_velocity','cell_volume') ) )
        px.append(reg.quantities['WeightedAverageQuantity']('momentum_x','cell_volume'))
        py.append(reg.quantities['WeightedAverageQuantity']('momentum_y','cell_volume'))
        pz.append(reg.quantities['WeightedAverageQuantity']('momentum_z','cell_volume'))
        ex.append(reg.quantities['WeightedAverageQuantity']('eng_x','cell_volume'))
        ey.append(reg.quantities['WeightedAverageQuantity']('eng_y','cell_volume'))
        ez.append(reg.quantities['WeightedAverageQuantity']('eng_z','cell_volume'))
    plt.clf()
    plt.plot(t,ex,label='ex')
    plt.plot(t,ey,label='ey')
    plt.plot(t,ez,label='ez')
    plt.legend(loc=0)
    plt.savefig('p42_Quan_%s_eng.pdf'%car.name)

    plt.clf()
    plt.plot(t,px,label='px')
    plt.plot(t,py,label='py')
    plt.plot(t,pz,label='pz')
    plt.legend(loc=0)
    plt.savefig('p42_Quan_%s_mom.pdf'%car.name)

    plt.clf()
    plt.plot(t,vx,label='vx')
    plt.plot(t,vy,label='vy')
    plt.plot(t,vz,label='vz')
    plt.plot(t,mach,label='mach')
    plt.legend(loc=0)
    plt.savefig('p42_Quan_%s_vel.pdf'%car.name)
