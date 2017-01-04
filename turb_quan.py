import taxi


#fleet = [aw11]
class quan_box():
    def __init__(self,car):
        self.name = car.name
        self.keys=['ex','ey','ez','vx','vy','vz','px','py','pz','mach']
        self.magkeys=['Bx','By','Bz','bx','by','bz','bx2','by2','bz2','beta','AlfMach','AlfvenSpeed']
    def __getitem__(self,key):
        return self.__dict__[key]
    def show(self,format="%10.2e"):
        N = len(self.frames)
        keys_to_print = self.keys
        if hasattr(self,'Bx'):
            keys_to_print += self.magkeys
        for key in keys_to_print:
            print "%5s"%key,
            print format*N%tuple(self[key])



    def __call__(self,car, tdyn=1, frames=None, extant_quan=None):
        #car.plot()
        print car.name
        if hasattr(extant_quan,'vx'):
            self.vx=extant_quan.vx
        else:
            self.vx=[]
        if hasattr(extant_quan,'vy'):
            self.vy=extant_quan.vy
        else:
            self.vy=[]
        if hasattr(extant_quan,'vz'):
            self.vz=extant_quan.vz
        else:
            self.vz=[]
        if hasattr(extant_quan,'mach'):
            self.mach=extant_quan.mach
        else:
            self.mach=[]
        if hasattr(extant_quan,'px'):
            self.px=extant_quan.px
        else:
            self.px=[]
        if hasattr(extant_quan,'py'):
            self.py=extant_quan.py
        else:
            self.py=[]
        if hasattr(extant_quan,'pz'):
            self.pz=extant_quan.pz
        else:
            self.pz=[]
        if hasattr(extant_quan,'ex'):
            self.ex=extant_quan.ex
        else:
            self.ex=[]
        if hasattr(extant_quan,'ey'):
            self.ey=extant_quan.ey
        else:
            self.ey=[]
        if hasattr(extant_quan,'ez'):
            self.ez=extant_quan.ez
        else:
            self.ez=[]
        if hasattr(extant_quan,'t'):
            self.t=extant_quan.t 
        else:
            self.t=[]

        if hasattr(extant_quan,'bx'):
            self.bx=extant_quan.bx
        else:
            self.bx=[]
        if hasattr(extant_quan,'by'):
            self.by=extant_quan.by
        else:
            self.by=[]
        if hasattr(extant_quan,'bz'):
            self.bz=extant_quan.bz
        else:
            self.bz=[]

        if hasattr(extant_quan,'bx2'):
            self.bx2=extant_quan.bx2
        else:
            self.bx2=[]
        if hasattr(extant_quan,'by2'):
            self.by2=extant_quan.by2
        else:
            self.by2=[]
        if hasattr(extant_quan,'bz2'):
            self.bz2=extant_quan.bz2
        else:
            self.bz2=[]


        if hasattr(extant_quan,'Bx'):
            self.Bx=extant_quan.Bx
        else:
            self.Bx=[]
        if hasattr(extant_quan,'By'):
            self.By=extant_quan.By
        else:
            self.By=[]
        if hasattr(extant_quan,'Bz'):
            self.Bz=extant_quan.Bz
        else:
            self.Bz=[]

        if hasattr(extant_quan,'Bfield_strength'):
            self.Bfield_strength=extant_quan.Bfield_strength
        else:
            self.Bfield_strength=[]

        if hasattr(extant_quan,'AlfMach'):
            self.AlfMach=extant_quan.AlfMach 
        else:
            self.AlfMach=[]
        if hasattr(extant_quan,'beta'):
            self.beta=extant_quan.beta 
        else:
            self.beta=[]
        if hasattr(extant_quan,'AlfvenSpeed'):
            self.AlfvenSpeed=extant_quan.AlfvenSpeed 
        else:
            self.AlfvenSpeed=[]
        if hasattr(extant_quan,'frames'):
            self.frames=extant_quan.frames
        else:
            self.frames=[]
        if hasattr(extant_quan,'tdyn'):
            self.tdyn = extant_quan.tdyn
        else:
            self.tdyn=tdyn

        if frames is None:
            frames = car.frames

        car.fill(frames[0])
        for frame in frames:
            if frame in self.frames:
                continue
            print "QUAN ON ", car.name, frame
            self.frames.append(frame)
            car.fill(frame)
            ds=car.ds
            self.t.append(car.ds.current_time/(tdyn) )
            car.region_type='all'
            reg = car.get_region(frame)
            total_volume = reg['cell_volume'].sum()
            volume = reg['cell_volume']
            self.mach.append( np.sqrt(reg.quantities['WeightedAverageQuantity']('mean_square_velocity','cell_volume').v ) )
            del reg['mean_square_velocity']
            self.vx.append(reg.quantities['WeightedAverageQuantity']('velocity_x','cell_volume').v)
            self.vy.append(reg.quantities['WeightedAverageQuantity']('velocity_y','cell_volume').v)
            self.vz.append(reg.quantities['WeightedAverageQuantity']('velocity_z','cell_volume').v)
            self.px.append(reg.quantities['WeightedAverageQuantity']('momentum_x','cell_volume').v)
            self.py.append(reg.quantities['WeightedAverageQuantity']('momentum_y','cell_volume').v)
            self.pz.append(reg.quantities['WeightedAverageQuantity']('momentum_z','cell_volume').v)
            for field in ['momentum_x','momentum_y','momentum_z']:
                del reg[field]
            self.ex.append(reg.quantities['WeightedAverageQuantity']('eng_x','cell_volume').v)
            self.ey.append(reg.quantities['WeightedAverageQuantity']('eng_y','cell_volume').v)
            self.ez.append(reg.quantities['WeightedAverageQuantity']('eng_z','cell_volume').v)

            if car.ds['HydroMethod'] in [4,6]:
                self.Bx.append( (reg['Bx']*volume).sum()/total_volume)
                self.By.append( (reg['By']*volume).sum()/total_volume)
                self.Bz.append( (reg['Bz']*volume).sum()/total_volume)

                self.bx.append( ( (reg['Bx']-self.Bx[-1])*volume).sum()/total_volume)
                self.by.append( ( (reg['By']-self.By[-1])*volume).sum()/total_volume)
                self.bz.append( ( (reg['Bz']-self.Bz[-1])*volume).sum()/total_volume)

                self.bx2.append( np.sqrt(( (reg['Bx']-self.Bx[-1])**2*volume).sum()/total_volume) )
                self.by2.append( np.sqrt(( (reg['By']-self.By[-1])**2*volume).sum()/total_volume) )
                self.bz2.append( np.sqrt(( (reg['Bz']-self.Bz[-1])**2*volume).sum()/total_volume) )

                self.Bfield_strength.append( (reg['magnetic_field_strength']*volume).sum()/total_volume)
                self.AlfvenSpeed.append( (volume*reg['magnetic_field_strength']/np.sqrt(np.pi*4*reg['density']) ).sum()/total_volume)

                self.AlfMach.append( self.mach[-1]/self.AlfvenSpeed[-1])
                self.beta.append( 2*(self.AlfMach[-1]/self.mach[-1])**2 )



        if car.ds['HydroMethod'] in [4,6]:
            plt.clf()
            plt.plot(self.t,self.Bx,label='Bx')
            plt.plot(self.t,self.By,label='By')
            plt.plot(self.t,self.Bz,label='Bz')
            plt.plot(self.t,self.bx,label='bx')
            plt.plot(self.t,self.by,label='by')
            plt.plot(self.t,self.bz,label='bz')
            plt.plot(self.t,self.bx2,label='bx rms')
            plt.plot(self.t,self.by2,label='by rms')
            plt.plot(self.t,self.bz2,label='bz rms')
            plt.legend(loc=0)
            plt.xlabel('t/tcode')
            plt.savefig('p42_Quan_%s_field_strength.pdf'%car.name)

            plt.clf()
            plt.plot(self.t, self.AlfMach, label='MA')
            plt.plot(self.t, self.mach, label="M")
            plt.plot(self.t, self.AlfvenSpeed,label='Va')
            plt.plot(self.t, self.beta, label='beta')
            plt.xlabel('t/tcode')
            plt.legend(loc=0)
            plt.savefig('p42_Quan_%s_MaM.pdf'%car.name)




        plt.clf()
        plt.plot(self.t,self.ex,label='ex')
        plt.plot(self.t,self.ey,label='ey')
        plt.plot(self.t,self.ez,label='ez')
        plt.legend(loc=0)
        plt.savefig('p42_Quan_%s_eng.pdf'%car.name)

        plt.clf()
        plt.plot(self.t,self.px,label='px')
        plt.plot(self.t,self.py,label='py')
        plt.plot(self.t,self.pz,label='pz')
        plt.legend(loc=0)
        plt.savefig('p42_Quan_%s_mom.pdf'%car.name)

        plt.clf()
        plt.plot(self.t,self.vx,label='vx')
        plt.plot(self.t,self.vy,label='vy')
        plt.plot(self.t,self.vz,label='vz')
        plt.plot(self.t,self.mach,label='mach')
        plt.legend(loc=0)
        plt.savefig('p42_Quan_%s_vel.pdf'%car.name)
