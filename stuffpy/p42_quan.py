import taxi


#fleet = [aw11]
class quan_box():
    def __init__(self,car):
        self.name = car.name
    def __call__(self,car, tdyn=1, frames=None):
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
        self.frames=[]

        
        if frames is None:
            frames = car.frames
        for frame in frames:
            self.frames.append(frame)
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
                beta.append( 2*(AlfMach[-1]/mach[-1])**2 )

            fptr = open("p42_aq42_quan.txt",'a')
            values  =    (t[-1],vx[-1],vy[-1],vz[-1],mach[-1],px[-1],py[-1],pz[-1],ex[-1],ey[-1],ez[-1]) 
            nfields = len(values)
            stringout =  "%0.16e "*nfields%tuple(values)
            fptr.write(stringout+"\n")
        self.Bfield_strength=  Bfield_strength
        self.vx=vx
        self.vy=vy
        self.vz=vz
        self.mach=mach
        self.px=px
        self.py=py
        self.pz=pz
        self.ex=ex
        self.ey=ey
        self.ez=ez
        self.t =t 
                                                 
        self.bx=bx
        self.by=by
        self.bz=bz
                                               
        self.bx2=bx2
        self.by2=by2
        self.bz2=bz2
                                                 
                                                 
        self.Bx=Bx
        self.By=By
        self.Bz=Bz
                                                 
                                                 
        self.AlfMach =AlfMach 
        self.beta =beta 
        self.AlfvenSpeed =AlfvenSpeed 

        self.tdyn = tdyn

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
            plt.xlabel('t/tcode')
            plt.savefig('p42_Quan_%s_field_strength.pdf'%car.name)

            plt.clf()
            plt.plot(t, AlfMach, label='MA')
            plt.plot(t, mach, label="M")
            plt.plot(t, AlfvenSpeed,label='Va')
            plt.plot(t, beta, label='beta')
            plt.xlabel('t/tcode')
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
