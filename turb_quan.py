import taxi
import p49_fields
import astropy.io.fits as pyfits
import glob
import fPickle
import os

#fleet = [aw11]
class quan_box():
    def __init__(self,car, plot_format='png'):
        self.car = car
        self.name = car.name
        self.keys=['ex','ey','ez','vx','vy','vz','px','py','pz','mach']
        self.magkeys=['Bx','By','Bz','bx','by','bz','bx2','by2','bz2','beta','AlfMach','AlfvenSpeed']
        self.all_fields = ['vx','vy','vz','mach','px','py','pz','ex','ey','ez','t','bx','by','bz','bx2','by2','bz2']
        self.all_fields +=['Bx','By','Bz','Bfield_strength','AlfMach','beta','AlfvenSpeed','frames']
        self.all_fields +=['ke_tot','ke_rel','grav_pot','grav_pot_2','gas_work']
        self.potential_written = False
        self.clobber=False
        self.stuff={}
        self.plot_format=plot_format
        #self.EBSlopePower={}

    def dump(self, pickle_name=None):
        if pickle_name is None:
            pickle_name = 'quan_box_%s.pickle'%self.car.name
        fPickle.dump(self.stuff,pickle_name)

    def load(self, pickle_name=None):
        if pickle_name is None:
            pickle_name = 'quan_box_%s.pickle'%self.car.name
        if len(glob.glob(pickle_name)):
            self.stuff = fPickle.load(pickle_name)
        else:
            print "NO SUCH FILE", pickle_name

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



    def plot_eb_vs_stuff(self,eb_quantity, axis, other_quantity):
        quad_frames = self.stuff['frames']
        eb_frames = self.stuff['EB'].keys()
        all_frames = np.unique(nar(quad_frames + eb_frames))


    def EBall(self,frames=None):
        if frames is None:
            frames = self.car.return_frames()
        for frame in frames:
            self.make_frbs(frame)
            self.QUEB(frame)
            self.EBslopes(frame)
        self.dump()
    def EBslopes(self,frame):
        import p49_slopes_powers
        reload(p49_slopes_powers)
        EBSlopePower=p49_slopes_powers.slopes_powers(self.car,frame, plot_format=self.plot_format)
        if not self.stuff.has_key('EB'):
            self.stuff['EB']={}
        self.stuff['EB'][frame]=EBSlopePower
    def QUEB(self, frame):
        import p49_QU2EB
        reload (p49_QU2EB)
        ds = self.car.load(frame)
        frb_dir = "%s/FRBs/"%self.car.directory
        p49_QU2EB.QU2EB(frb_dir,frame)
#        if frames is None:
#            frames = self.car.return_frames()
#        for frame in frames:
#            self.make_frbs(frame)
#            self.fit_slopes()
#            self.plot_eebb()
    def make_frbs(self,frame):
        fields=[]
        for axis in ['x', 'y', 'z']:
          n0=1; p=1 #n0 in [19,39,1945] and p=0
          fields.append( (axis,'Q%s_n0-%04d_p-%d'%(axis,n0,p))   )
          fields.append( (axis,'U%s_n0-%04d_p-%d'%(axis,n0,p))   )
        ds = self.car.load(frame)
        res = ds.parameters['TopGridDimensions'][2 + ord('x') - ord(axis)] # zyx order

        for axis, field in fields :
            outputdir = "%s/FRBs/"%self.car.directory
            if not os.access(outputdir, os.F_OK):
                os.mkdir(outputdir)
            outfile = outputdir+"/DD%.4d_%s.fits" %(frame,field)
            if os.access(outfile, os.F_OK) and not self.clobber:
                print "FRB exists: %s"%outfile
            else:
                print "FRB being produced: %s"%outfile
                res = ds.parameters['TopGridDimensions'][2 + ord('x') - ord(axis)]
                proj = ds.proj(field,axis)
                frb = proj.to_frb(1,res)
                hdu = pyfits.PrimaryHDU(frb[field])
                hdulist = pyfits.HDUList([hdu])
                hdulist.writeto(outfile,clobber=True)
                print "wrote", outfile



        #car.plot()
    def __call__(self, tdyn=1, frames=None):
        if frames is None:
            frames = self.car.return_frames()
        for frame in self.car.return_frames():
            self.quadratic(self.car,frame,tdyn)
        #car.plot()
    def quadratic(self,car, frame, tdyn=1):
        print car.name

        for field in self.all_fields:
            if not self.stuff.has_key(field):
                self.stuff[field] = []

        if self.stuff.has_key('tdyn'):
            self.tdyn = self.stuff['tdyn']
        else:
            self.tdyn=tdyn
        self.stuff['tdyn']=self.tdyn



        car.region_type='all'
        ds=car.load(frame)
        self.potential_written = False
        if car.ds['SelfGravity']:
            car.ds.create_field_info()
            self.potential_written = 'PotentialField' in [k[1] for k in car.ds.field_info.keys()]

        if frame in self.stuff['frames'] and self.clobber == False:
            return
        print "QUAN ON ", car.name, frame
        self.stuff['frames'].append(frame)
        self.stuff['t'].append(car.ds['InitialTime']/(tdyn) )
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
        if self.potential_written:
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
        #picklename = "turb_quan_temp_%s.pickle"%car.name
        self.dump()
        #self.plot()



    def plot(self):
        print "here's the thing"
        car = self.car
        ds=car.load()
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
            outname = '%s_quan_field_strength.%s'%(car.name,self.plot_format)
            plt.savefig(outname)
            print outname

            plt.clf()
            plt.plot(self.stuff['t'], self.stuff['AlfMach'], label='MA')
            plt.plot(self.stuff['t'], self.stuff['mach'], label="M")
            plt.plot(self.stuff['t'], self.stuff['AlfvenSpeed'],label='Va')
            plt.plot(self.stuff['t'], self.stuff['beta'], label='beta')
            plt.xlabel('t/tcode')
            plt.legend(loc=0)
            outname = '%s_quan_MaM.%s'%(car.name,self.plot_format)
            plt.savefig(outname)
            print outname
            plt.yscale('log')
            outname = '%s_quan_MaM_log.%s'%(car.name,self.plot_format)
            plt.savefig(outname)
            print outname




        plt.clf()
        plt.plot(self.stuff['t'],self.stuff['ex'],label='ex')
        plt.plot(self.stuff['t'],self.stuff['ey'],label='ey')
        plt.plot(self.stuff['t'],self.stuff['ez'],label='ez')
        plt.legend(loc=0)
        outname = '%s_quan_eng.%s'%(car.name,self.plot_format)
        plt.savefig(outname)
        print outname

        plt.clf()
        plt.plot(self.stuff['t'],self.stuff['px'],label='px')
        plt.plot(self.stuff['t'],self.stuff['py'],label='py')
        plt.plot(self.stuff['t'],self.stuff['pz'],label='pz')
        plt.legend(loc=0)
        outname = '%s_quan_mom.%s'%(car.name,self.plot_format)
        plt.savefig(outname)
        print outname

        plt.clf()
        plt.plot(self.stuff['t'],self.stuff['vx'],label='vx')
        plt.plot(self.stuff['t'],self.stuff['vy'],label='vy')
        plt.plot(self.stuff['t'],self.stuff['vz'],label='vz')
        plt.plot(self.stuff['t'],self.stuff['mach'],label='mach')
        plt.legend(loc=0)
        outname = '%s_quan_vel.%s'%(car.name,self.plot_format)
        plt.savefig(outname)
        print outname

        plt.clf()
        plt.plot(self.stuff['t'],self.stuff['ke_tot'],label='ke_tot')
        plt.plot(self.stuff['t'],self.stuff['ke_rel'],label='ke_rel')
        if self.potential_written: 
            plt.plot(self.stuff['t'],self.stuff['grav_pot'],label='grav_pot')
        plt.plot(self.stuff['t'],self.stuff['gas_work'],label='gas_work')
        plt.legend(loc=0)
        outname = '%s_quan_energy.%s'%(car.name,self.plot_format)
        plt.savefig(outname)
        print outname
