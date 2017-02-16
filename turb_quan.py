import taxi


#fleet = [aw11]
class quan_box():
    def __init__(self,car):
        self.name = car.name
        self.keys=['ex','ey','ez','vx','vy','vz','px','py','pz','mach']
        self.magkeys=['Bx','By','Bz','bx','by','bz','bx2','by2','bz2','beta','AlfMach','AlfvenSpeed']
        self.stuff={}
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

        all_fields = ['vx','vy','vz','mach','px','py','pz','ex','ey','ez','t','bx','by','bz','bx2','by2','bz2']
        all_fields +=['Bx','By','Bz','Bfield_strength','AlfMach','beta','AlfvenSpeed','frames']
        all_fields +=['ke_tot','ke_rel','grav_pot','grav_pot_2','gas_work']
        for field in all_fields:
            if extant_quan is not None and extant_quan.stuff.has_key(field):
                self.stuff[field] = extant_quan[field]
            else:
                self.stuff[field] = []

        if hasattr(extant_quan,'tdyn'):
            self.tdyn = extant_quan.tdyn
        else:
            self.tdyn=tdyn

        if frames is None:
            frames = car.frames

        car.fill(frames[0])
        potential_written = False
        if car.ds['SelfGravity']:
            car.ds.create_field_info()
            potential_written = 'PotentialField' in [k[1] for k in car.ds.field_info.keys()]

        for frame in frames:
            if frame in self.stuff['frames']:
                continue
            print "QUAN ON ", car.name, frame
            self.stuff['frames'].append(frame)
            car.fill(frame)
            ds=car.ds
            self.stuff['t'].append(car.ds['InitialTime']/(tdyn) )
            car.region_type='all'
            reg = car.get_region(frame)
            total_volume = reg['cell_volume'].sum()
            volume = reg['cell_volume']
            self.stuff['mach'].append( np.sqrt(reg.quantities['WeightedAverageQuantity']('mean_square_velocity','cell_volume').in_units('code_velocity**2').v ) )
            self.stuff['vx'].append(reg.quantities['WeightedAverageQuantity']('velocity_x','cell_volume').in_units('code_velocity').v)
            self.stuff['vy'].append(reg.quantities['WeightedAverageQuantity']('velocity_y','cell_volume').in_units('code_velocity').v)
            self.stuff['vz'].append(reg.quantities['WeightedAverageQuantity']('velocity_z','cell_volume').in_units('code_velocity').v)
            self.stuff['px'].append(reg.quantities['WeightedAverageQuantity']('momentum_x','cell_volume').in_units('code_density*code_velocity').v)
            self.stuff['py'].append(reg.quantities['WeightedAverageQuantity']('momentum_y','cell_volume').in_units('code_density*code_velocity').v)
            self.stuff['pz'].append(reg.quantities['WeightedAverageQuantity']('momentum_z','cell_volume').in_units('code_density*code_velocity').v)
            self.stuff['ex'].append(reg.quantities['WeightedAverageQuantity']('eng_x','cell_volume').in_units('code_density*code_velocity**2').v)
            self.stuff['ey'].append(reg.quantities['WeightedAverageQuantity']('eng_y','cell_volume').in_units('code_density*code_velocity**2').v)
            self.stuff['ez'].append(reg.quantities['WeightedAverageQuantity']('eng_z','cell_volume').in_units('code_density*code_velocity**2').v)
            self.stuff['ke_tot'].append(reg.quantities['WeightedAverageQuantity']('kinetic_energy','cell_volume').in_units('code_density*code_velocity**2').v)
            reg.set_field_parameter('bulk_velocity',ds.arr([self.stuff['vx'][-1],self.stuff['vy'][-1], self.stuff['vz'][-1]],'code_velocity'))
            self.stuff['ke_rel'].append(reg.quantities['WeightedAverageQuantity']('rel_kinetic_energy','cell_volume').in_units('code_density*code_velocity**2').v)
            if potential_written:
                self.stuff['grav_pot'].append(reg.quantities['WeightedAverageQuantity']('grav_pot','cell_volume').in_units('code_density*code_velocity**2').v)
            self.stuff['gas_work'].append(reg.quantities['WeightedAverageQuantity']('gas_work','cell_volume').in_units('code_density*code_velocity**2').v)

            if car.ds['HydroMethod'] in [4,6]:
                self.stuff['Bx'].append( (reg['Bx']*volume).sum()/total_volume)
                self.stuff['By'].append( (reg['By']*volume).sum()/total_volume)
                self.stuff['Bz'].append( (reg['Bz']*volume).sum()/total_volume)

                self.stuff['bx'].append( ( (reg['Bx']-self.stuff['Bx'][-1])*volume).sum()/total_volume)
                self.stuff['by'].append( ( (reg['By']-self.stuff['By'][-1])*volume).sum()/total_volume)
                self.stuff['bz'].append( ( (reg['Bz']-self.stuff['Bz'][-1])*volume).sum()/total_volume)

                self.stuff['bx2'].append( np.sqrt(( (reg['Bx']-self.stuff['Bx'][-1])**2*volume).sum()/total_volume) )
                self.stuff['by2'].append( np.sqrt(( (reg['By']-self.stuff['By'][-1])**2*volume).sum()/total_volume) )
                self.stuff['bz2'].append( np.sqrt(( (reg['Bz']-self.stuff['Bz'][-1])**2*volume).sum()/total_volume) )

                self.stuff['Bfield_strength'].append( (reg['magnetic_field_strength']*volume).sum()/total_volume)
                self.stuff['AlfvenSpeed'].append( (volume*reg['magnetic_field_strength']/np.sqrt(np.pi*4*reg['density']) ).sum()/total_volume)

                self.stuff['AlfMach'].append( self.stuff['mach'][-1]/self.stuff['AlfvenSpeed'][-1])
                self.stuff['beta'].append( 2*(self.stuff['AlfMach'][-1]/self.stuff['mach'][-1])**2 )



        if car.ds['HydroMethod'] in [4,6]:
            plt.clf()
            plt.plot(self.stuff['t'],self.stuff['Bx'],label='Bx')
            plt.plot(self.stuff['t'],self.stuff['By'],label='By')
            plt.plot(self.stuff['t'],self.stuff['Bz'],label='Bz')
            plt.plot(self.stuff['t'],self.stuff['bx'],label='bx')
            plt.plot(self.stuff['t'],self.stuff['by'],label='by')
            plt.plot(self.stuff['t'],self.stuff['bz'],label='bz')
            plt.plot(self.stuff['t'],self.stuff['bx2'],label='bx rms')
            plt.plot(self.stuff['t'],self.stuff['by2'],label='by rms')
            plt.plot(self.stuff['t'],self.stuff['bz2'],label='bz rms')
            plt.legend(loc=0)
            plt.xlabel('t/tcode')
            plt.savefig('p42_Quan_%s_field_strength.pdf'%car.name)

            plt.clf()
            plt.plot(self.stuff['t'], self.stuff['AlfMach'], label='MA')
            plt.plot(self.stuff['t'], self.stuff['mach'], label="M")
            plt.plot(self.stuff['t'], self.stuff['AlfvenSpeed'],label='Va')
            plt.plot(self.stuff['t'], self.stuff['beta'], label='beta')
            plt.xlabel('t/tcode')
            plt.legend(loc=0)
            plt.savefig('p42_Quan_%s_MaM.pdf'%car.name)




        plt.clf()
        plt.plot(self.stuff['t'],self.stuff['ex'],label='ex')
        plt.plot(self.stuff['t'],self.stuff['ey'],label='ey')
        plt.plot(self.stuff['t'],self.stuff['ez'],label='ez')
        plt.legend(loc=0)
        plt.savefig('p42_Quan_%s_eng.pdf'%car.name)

        plt.clf()
        plt.plot(self.stuff['t'],self.stuff['px'],label='px')
        plt.plot(self.stuff['t'],self.stuff['py'],label='py')
        plt.plot(self.stuff['t'],self.stuff['pz'],label='pz')
        plt.legend(loc=0)
        plt.savefig('p42_Quan_%s_mom.pdf'%car.name)

        plt.clf()
        plt.plot(self.stuff['t'],self.stuff['vx'],label='vx')
        plt.plot(self.stuff['t'],self.stuff['vy'],label='vy')
        plt.plot(self.stuff['t'],self.stuff['vz'],label='vz')
        plt.plot(self.stuff['t'],self.stuff['mach'],label='mach')
        plt.legend(loc=0)
        plt.savefig('p42_Quan_%s_vel.pdf'%car.name)

        plt.clf()
        plt.plot(self.stuff['t'],self.stuff['ke_tot'],label='ke_tot')
        plt.plot(self.stuff['t'],self.stuff['ke_rel'],label='ke_rel')
        if potential_written: 
            plt.plot(self.stuff['t'],self.stuff['grav_pot'],label='grav_pot')
        plt.plot(self.stuff['t'],self.stuff['gas_work'],label='gas_work')
        plt.legend(loc=0)
        plt.savefig('p42_Quan_%s_energy.pdf'%car.name)
