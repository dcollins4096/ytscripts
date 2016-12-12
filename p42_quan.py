import taxi
aq32 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq32_m2.9_drive1_32',name='aq32',frames=range(0,42,2),fields=['density'])
aq31 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq31_m2.9_drive0.5_32',name='aq31',frames=range(0,42,2),fields=['density'])
aq16 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq16_ppm_m9_drive0_noamr_128',name='aq15',frames=range(17),fields=['density'])
aq33 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq33_m9_drive0.5_32',name='aq33',frames=range(0,110,10))
aq34 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq34_m9_drive0_p60_32', name='aq34', frames=range(0,110,10))
aq38 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq38_m2.9_drive0_32_p59_ppm',name='aq38',frames=range(0,42,2),fields=['density'])
aq39 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq39_m2.9_drive0_32_p59_ppm' ,name='aq39',frames=range(0,42,2),fields=['density'])
aq40 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq40_m2.9_drive0_32_p59_hydro6_newAddVel' ,name='aq40',frames=range(0,42,2),fields=['density'])
aq41 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq41_m9_drive0.5_32_p59_ppm' ,name='aq41',frames=range(0,42,2),fields=['density'])
aq42 = taxi.taxi(directory='/scratch/00369/tg456484/Paper42_NewAK/aq42_m9_drive0.5_512_p59_ppm' ,name='aq42',frames=range(90,0,-10),fields=['density'])
aw11 = taxi.taxi(directory='/scratch/00369/tg456484/Paper49_EB/aw11_M1_MA1_256_p59' ,name='aw11',frames=range(0,210,10),fields=['density'])


fleet = [aw11]
for car in fleet:
    #car.plot()
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

    bx=[]
    by=[]
    bz=[]

    bx2=[]
    by2=[]
    bz2=[]


    Bx=[]
    By=[]
    Bz=[]

    Bfield_strength=[]

    AlfMach = []
    beta = []
    AlfvenSpeed = []

    tdyn = 0.5/9
    
    for frame in car.frames:
        car.fill(frame)
        ds=car.ds
        t.append(car.ds.current_time/(tdyn) )
        car.region_type='all'
        reg = car.get_region(frame)
        total_volume = reg['cell_volume'].sum()
        volume = reg['cell_volume']
        mach.append( np.sqrt(reg.quantities['WeightedAverageQuantity']('mean_square_velocity','cell_volume').v ) )
        vx.append(reg.quantities['WeightedAverageQuantity']('velocity_x','cell_volume').v)
        vy.append(reg.quantities['WeightedAverageQuantity']('velocity_y','cell_volume').v)
        vz.append(reg.quantities['WeightedAverageQuantity']('velocity_z','cell_volume').v)
        px.append(reg.quantities['WeightedAverageQuantity']('momentum_x','cell_volume').v)
        py.append(reg.quantities['WeightedAverageQuantity']('momentum_y','cell_volume').v)
        pz.append(reg.quantities['WeightedAverageQuantity']('momentum_z','cell_volume').v)
        ex.append(reg.quantities['WeightedAverageQuantity']('eng_x','cell_volume').v)
        ey.append(reg.quantities['WeightedAverageQuantity']('eng_y','cell_volume').v)
        ez.append(reg.quantities['WeightedAverageQuantity']('eng_z','cell_volume').v)

        if car.ds['HydroMethod'] in [4,6]:
            Bx.append( (reg['Bx']*volume).sum()/total_volume)
            By.append( (reg['By']*volume).sum()/total_volume)
            Bz.append( (reg['Bz']*volume).sum()/total_volume)

            bx.append( ( (reg['Bx']-Bx[-1])*volume).sum()/total_volume)
            by.append( ( (reg['By']-By[-1])*volume).sum()/total_volume)
            bz.append( ( (reg['Bz']-Bz[-1])*volume).sum()/total_volume)

            bx2.append( np.sqrt(( (reg['Bx']-Bx[-1])**2*volume).sum()/total_volume) )
            by2.append( np.sqrt(( (reg['By']-By[-1])**2*volume).sum()/total_volume) )
            bz2.append( np.sqrt(( (reg['Bz']-Bz[-1])**2*volume).sum()/total_volume) )

            Bfield_strength.append( (reg['magnetic_field_strength']*volume).sum()/total_volume)
            AlfvenSpeed.append( (volume*reg['magnetic_field_strength']/np.sqrt(np.pi*4*reg['density']) ).sum()/total_volume)

            AlfMach.append( mach[-1]/AlfvenSpeed[-1])
            beta.append( (mach[-1]/AlfMach[-1])**2 )

        fptr = open("p42_aq42_quan.txt",'a')
        values  =    (t[-1],vx[-1],vy[-1],vz[-1],mach[-1],px[-1],py[-1],pz[-1],ex[-1],ey[-1],ez[-1]) 
        nfields = len(values)
        stringout =  "%0.16e "*nfields%tuple(values)
        fptr.write(stringout+"\n")

    if car.ds['HydroMethod'] in [4,6]:
        plt.clf()
        plt.plot(t,Bx,label='Bx')
        plt.plot(t,By,label='By')
        plt.plot(t,Bz,label='Bz')
        plt.plot(t,bx,label='bx')
        plt.plot(t,by,label='by')
        plt.plot(t,bz,label='bz')
        plt.plot(t,bx2,label='bx rms')
        plt.plot(t,by2,label='by rms')
        plt.plot(t,bz2,label='bz rms')
        plt.legend(loc=0)
        plt.xlabel('t/tdyn')
        plt.savefig('p42_Quan_%s_field_strength.pdf'%car.name)

        plt.clf()
        plt.plot(t, AlfMach, label='MA')
        plt.plot(t, mach, label="M")
        plt.plot(t, AlfvenSpeed,label='Va')
        plt.plot(t, beta, label='beta')
        plt.xlabel('t/tdyn')
        plt.legend(loc=0)
        plt.savefig('p42_Quan_%s_MaM.pdf'%car.name)




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
