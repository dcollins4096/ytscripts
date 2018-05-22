execfile('go_lite')
import taxi
import p49_fields
import p49_labels
import p49_QU2EB
import xtra_energy_fields
import astropy.io.fits as pyfits
import glob
import fPickle
import os
import numpy as np
import time
import glob
import h5py
import re
nar = np.array


frbname = p49_QU2EB.frbname
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

    def merge(self,dict2):
        dict1 = self.stuff
        #dict2 = other_quan.stuff
        frames1 = dict1['frames']
        frames2 = dict2['frames']
        for i, frame in enumerate(frames2):
            if frame not in frames1:
                for key in dict1:
                    if key in [ 'EB', 'grav_pot_2', 'tdyn']:
                        continue
                    if key in ['grav_pot'] and not self.potential_written:
                        continue
                    dict1[key].append(dict2[key][i])

        if 'EB' in dict2:
            if 'EB' not in dict1:
                dict1['EB'] = {}
            for frame in dict2['EB']:
                if frame not in dict1['EB']:
                    dict1['EB'][frame] = dict2['EB'][frame]
    def dump(self, pickle_name=None):
        #pdb.set_trace()
        if pickle_name is None:
            pickle_name = 'quan_box_%s.pickle'%self.car.name
        lock_name = pickle_name + ".lock"
        counter = 0
        while os.path.exists(lock_name) and counter < 5:
            print lock_name, "exists"
            counter += 1
            time.sleep(1)
        if counter > 15:
            pickle_name = pickle_name+"%d"%len(glob.glob("%s*"%pickle_name))
        fptr = open(lock_name,"w+")
        fptr.write("in use\n");
        fptr.close()
        if os.path.exists(pickle_name):
            other_pickle = fPickle.load(pickle_name)
            self.merge(other_pickle)
            
        fPickle.dump(self.stuff,pickle_name)
        os.remove(lock_name)


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
        if not self.stuff.has_key('EB'):
            self.stuff['EB']={}
        if frames is None:
            frames = self.car.return_frames()
        for frame in frames:
            if frame not in self.stuff['EB']:
                self.make_frbs(frame)
                self.QUEB(frame)
                self.EBslopes(frame)
        self.dump()
    def EBslopes(self,frame):
        EBSlopePower=p49_QU2EB.slopes_powers(self.car,frame, plot_format=self.plot_format)
        if not self.stuff.has_key('EB'):
            self.stuff['EB']={}
        self.stuff['EB'][frame]=EBSlopePower
    def GetQUEB(self,frame):

        frb_dir = "%s/"+frbname"+/"%self.car.directory
        Qlist = glob.glob(frb_dir+'/DD%04d_Q[xyz]*.fits'%frame)
        Ulist = []
        # write the output near the input
        self.QUEBarr = {'Q':{}, 'U':{}, 'E':{}, 'B':{}}
        for Qfile in Qlist:
            mo = re.match('(.*/DD[0-9]{4}_)Q([xyz].*)(.fits)',Qfile)
            Ufile = mo.group(1)+'U'+mo.group(2)+'.fits'
            Ulist.append(Ufile)
            mo = re.match('(.*/DD[0-9]{4}_)Q([xyz].*)(.fits)',Qfile)
            outroot = mo.group(1)
            outsuf = mo.group(2)
            Efile = outroot+'E'+outsuf+'.fits'
            Bfile = outroot+'B'+outsuf+'.fits'
            Clfile = outroot+'Cl'+outsuf+'.dat'
            self.QUEBarr['Q'][Qfile] = np.array(pyfits.open(Qfile)[0].data,dtype=np.double)
            self.QUEBarr['U'][Ufile] = np.array(pyfits.open(Ufile)[0].data,dtype=np.double)
            self.QUEBarr['E'][Efile] = np.array(pyfits.open(Efile)[0].data,dtype=np.double)
            self.QUEBarr['B'][Bfile] = np.array(pyfits.open(Bfile)[0].data,dtype=np.double)


    def QUEB(self, frame):
        #ds = self.car.load(frame)
        frb_dir = "%s/"+frbname"+/"%self.car.directory
        p49_QU2EB.QU2EB(frb_dir,frame)
#        if frames is None:
#            frames = self.car.return_frames()
#        for frame in frames:
#            self.make_frbs(frame)
#            self.fit_slopes()
#            self.plot_eebb()
    def make_frbs(self,frame, axes=['x','y','z']):
        fields=[]
        for axis in axes:
          n0=1; p=1 #n0 in [19,39,1945] and p=0
          #fields.append( (axis,'Q%s_n0-%04d_p-%d'%(axis,n0,p))   )
          #fields.append( (axis,'U%s_n0-%04d_p-%d'%(axis,n0,p))   )
          fields.append( (axis,'Q%s'%(axis))   )
          fields.append( (axis,'U%s'%(axis))   )
          fields.append( (axis,'density') )

        ds = None  #this is somewhat awkward, but useful for avoiding simulations
                   #  that have only products, not datasest
        for axis, field in fields :
            outputdir = "%s/"+frbname"+/"%self.car.directory
            if not os.access(outputdir, os.F_OK):
                os.mkdir(outputdir)
            #Hm.  Q and U have the name in the field, but others don't.
            if field[0] in 'QU' and field[1] in 'xyz':
                field_name = field
            else:
                field_name = field + "_"+axis
            outfile = outputdir+"/DD%.4d_%s.fits" %(frame,field_name)
            #move this in the 'make' conditional?
            #if ds is None:
            #    ds = self.car.load(frame)
            #    res = ds.parameters['TopGridDimensions'][2 + ord('x') - ord(axis)] # zyx order
            if os.access(outfile, os.F_OK) and not self.clobber:
                print "FRB exists: %s"%outfile
            else:
                if ds is None:
                    ds = self.car.load(frame)
                    res = ds.parameters['TopGridDimensions'][2 + ord('x') - ord(axis)] # zyx order
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

        if frame in self.stuff['frames'] and self.clobber == False:
            return
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

        print "QUAN ON ", car.name, frame
        self.stuff['frames'].append(frame)
        self.stuff['t'].append(car.ds['InitialTime']/(tdyn) )
        reg = car.get_region(frame)
        total_volume = reg['cell_volume'].sum()
        volume = reg['cell_volume']
        #self.stuff['mass'].append( (reg.quantities['WeightedAverageQuantity']('density', 'cell_volume')*total_volume).in_units('code_mass'))
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





    def plot(self, HydroMethod = None):
        print "Hydro Method", HydroMethod
        def nom_string(val):
            if val > 10 or val < 0.1:
                return "%0.1e"%val
            else:
                return "%0.1f"%val
        print "here's the thing"

        car = self.car
        if HydroMethod is None:
            ds=car.load()
            HydroMethod = car.ds['HydroMethod']
        tn = p49_labels.nominal.get(car.name,None)
        times = nar(self.stuff['t'])
        ar = np.argsort(times)
        times=times[ar]
        for thing in self.stuff:
            if thing in ['tdyn', 'EB']:
                continue
            if len(self.stuff[thing]) == 0:
                continue
            self.stuff[thing] = nar(self.stuff[thing])[ar]
        time_label  =r'$t [\rm{code}]$'
        if tn:
            time_label += r' $(t_{\rm{cross}}= %s)$'% nom_string(0.5/tn['mach'])

        sqrtfourpi=np.sqrt(4*np.pi)

        if HydroMethod in [4,6]:
            plt.clf()
            n_points = len(times)
            my_ones = np.ones(n_points)
            
            if tn is not None:
                nom_val = tn['field_cgs']/sqrtfourpi
                st = "nominal "+ (nom_string(nom_val))
                plt.plot(times, my_ones*nom_val,label=st,c=[0.5]*4)
            plt.plot(times,self.stuff['Bx'],label='Bx')
            plt.plot(times,self.stuff['By'],label='By')
            plt.plot(times,self.stuff['Bz'],label='Bz')
            plt.plot(times,self.stuff['bx'],label='bx')
            plt.plot(times,self.stuff['by'],label='by')
            plt.plot(times,self.stuff['bz'],label='bz')
            plt.plot(times,self.stuff['bx2'],label='bx rms')
            plt.plot(times,self.stuff['by2'],label='by rms')
            plt.plot(times,self.stuff['bz2'],label='bz rms')
            plt.legend(loc=0)
            plt.ylabel('MagneticField'); plt.xlabel(time_label)
            outname = '%s_quan_field_strength.%s'%(car.name,self.plot_format)
            plt.savefig(outname)
            print outname

            plt.clf()
            c='r'
            if tn is not None:
                plt.plot(times, my_ones*tn['AlfMach'],label='nominal',c=c,linestyle='--')
            plt.plot(times, self.stuff['AlfMach'], label='MA',c=c )
            c='g'
            if tn is not None:
                plt.plot(times, my_ones*tn['mach'],c=c,linestyle='--')
            plt.plot(times, self.stuff['mach'], label="M", c='g')
            c='c'
            if tn is not None:
                plt.plot(times, my_ones*(tn['mach']/tn['AlfMach']),c=c,linestyle='--')
            plt.plot(times, self.stuff['AlfvenSpeed'],label='Va',c=c)
            c='b'
            if tn is not None:
                plt.plot(times, my_ones*(10**tn['logbeta']),c=c,linestyle='--')
            plt.plot(times, self.stuff['beta'], label='beta',c=c)
            plt.ylabel('Dimensionless'); plt.xlabel(time_label)
            plt.legend(loc=0)
            outname = '%s_quan_MaM.%s'%(car.name,self.plot_format)
            plt.savefig(outname)
            print outname
            plt.yscale('log')
            outname = '%s_quan_MaM_log.%s'%(car.name,self.plot_format)
            plt.savefig(outname)
            print outname

        plt.clf()
        plt.plot(times,self.stuff['ex'],label='ex')
        plt.plot(times,self.stuff['ey'],label='ey')
        plt.plot(times,self.stuff['ez'],label='ez')
        plt.ylabel('Partial Energies'); plt.xlabel(time_label)
        plt.legend(loc=0)
        outname = '%s_quan_eng.%s'%(car.name,self.plot_format)
        plt.savefig(outname)
        print outname

        plt.clf()
        plt.plot(times,self.stuff['px'],label='px')
        plt.plot(times,self.stuff['py'],label='py')
        plt.plot(times,self.stuff['pz'],label='pz')
        plt.legend(loc=0)
        outname = '%s_quan_mom.%s'%(car.name,self.plot_format)
        plt.ylabel('Momentum'); plt.xlabel(time_label)
        plt.savefig(outname)
        print outname

        plt.clf()
        plt.plot(times,self.stuff['vx'],label='vx',c='r')
        plt.plot(times,self.stuff['vy'],label='vy',c='g')
        plt.plot(times,self.stuff['vz'],label='vz',c='b')
        c='k'
        if tn is not None:
            plt.plot(times, my_ones*tn['mach'],c=c,linestyle='--')
        plt.plot(times,self.stuff['mach'],label='mach',c=c)
        plt.ylabel('velocities'); plt.xlabel(time_label)
        plt.legend(loc=0)
        outname = '%s_quan_vel.%s'%(car.name,self.plot_format)
        plt.savefig(outname)
        print outname

        plt.clf()
        plt.plot(times,self.stuff['ke_tot'],label='ke_tot')
        plt.plot(times,self.stuff['ke_rel'],label='ke_rel')
        if self.potential_written: 
            plt.plot(times,self.stuff['grav_pot'],label='grav_pot')
        plt.plot(times,self.stuff['gas_work'],label='gas_work')
        plt.ylabel('Energies'); plt.xlabel(time_label)
        plt.legend(loc=0)
        outname = '%s_quan_energy.%s'%(car.name,self.plot_format)
        plt.savefig(outname)
        print outname
