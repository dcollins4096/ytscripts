#from matplotlib.backends.backend_agg import FigureCanvasAgg
#from matplotlib import use; use('Agg') 
import yt
import pyfits
import glob
import pdb
import os
import matplotlib.pyplot as plt
from astrodendro import Dendrogram

ef = execfile
ef('fields_overhaul.py')

def stat(array,strin='', format='%0.16e'):
    """print basic statistics.  Min and max and size."""
    template = '['+format+','+format+'] %s %s'
    print template%(array.min(),array.max(),array.shape,strin)


class access_thing():
    def __init__(self, sim, res, frame, n0='', p='', machine='nazare'):
        """Sets location of enzo data set and fits location
        sim: the simulation
        res: resolution of the simulation
        frame: the frame number
        n0: critical density where power law changes for stokes parameters
        p: power law index for density above n0
        """
        self.sim=sim
        self.frame=frame
        self.n0 = n0     # Critical density
        self.p = p       # power law exponent
        restart_sim_template = "%s/RS%04d/restart%04d"

        self.sim_dir = {'High_512':'B02/512','Mid_512':'B2/512','Low_512':'B20/512',
                        'High_256':'B02/256','Mid_256':'B2/256','Low_256':'B20/256'}[sim]        

        if machine == 'nazare':
            self.enzo_location = '/scratch1/dcollins/Paper08/%s' %self.sim_dir
            self.fits_location = '/scratch2/cdb09f/enzotest/%s_test/%sres' %(sim, res)
            #self.plot_location = '/scratch2/cdb09f/enzotest/%s_res_plots'%res 
            self.plot_location = '/scratch2/cdb09f/enzotest/%s_res_plots/EB_power_spectra'%res             
            self.plot_location_P_frac             = '%s/P_frac_plots'%self.plot_location 
            self.plot_location_field_strength     = '%s/field_strength_plots'%self.plot_location 
            self.plot_location_viewing_angle      = '%s/viewing_angle'%self.plot_location   
            self.plot_location_time_dependence    = '%s/time_dependence_plots'%self.plot_location
            self.plot_location_threshold_density  = '%s/threshold_density_plots'%self.plot_location    
            self.plot_location_leaves             = '%s/leaves_plots'%self.plot_location    
            self.plot_location_P_theta_vs_B_theta = '%s/P_theta_vs_B_theta'%self.plot_location
            self.plot_location_time_dependence    = '%s/time_dependence'%self.plot_location
            self.plot_location_leaf_contours      = '%s/leaf_contours'%self.plot_location
            self.plot_location_mass_vs_smoothness = '%s/mass_vs_smoothness'%self.plot_location 
            self.plot_location_QU2EB              = '%s/field_strength_and_frame'%self.plot_location                                                                                                                                                                                                                                                                                                                                                                                                                                                   
            #self.data_location = '%s/data/frame_%04d'%(self.fits_location, self.frame)
            self.data_location = '%s/data/frame_%04d/QU2EB'%(self.fits_location, self.frame)
            self.dendrogram_location = self.fits_location

        elif machine == 'corey_local':
            self.enzo_location = '/Users/Corey/Desktop/stuff/FSU/Research/enzodata/%s' %self.sim_dir
            self.fits_location = '/Users/Corey/Desktop/stuff/FSU/Research/%s_test/%sres' %(sim, res)
            self.plot_location = '%s/plots/frame_%04d'%(self.fits_location, self.frame)
            self.plot_location_leaf_contours = '/scratch2/cdb09f/enzotest/%s_res_plots/leaf_contours'%res                                                                                                                                                                                                                                        
            self.data_location = '%s/data/frame_%04d'%(self.fits_location, self.frame)            
            self.dendrogram_location = self.fits_location

        self.enzo_dataset_template = restart_sim_template
        
        self.enzo_dataset = self.enzo_dataset_template%(self.enzo_location, frame, frame)
          

    def get_fits_name(self,field,axis):
        #self.fits_path = "%s/fits/frame_%04d/%s"%(self.fits_location, self.frame, axis)
        self.fits_path = "%s/fits/frame_%04d/removed_N_avg-NH_Av_from_epsilon/%s"%(self.fits_location, self.frame, axis)
        self.B_fits_path = "%s/fits/frame_%04d"%(self.fits_location, self.frame)
        self.density_fits_name = "%s/%s.fits"%(self.fits_path, field)

        if field == 'Bx' or field == 'By' or field == 'Bz':
            self.fits_name = "%s/%s.fits"%(self.B_fits_path, field)
        else:
            self.fits_name = "%s/%s.fits"%(self.fits_path, field)

    def get_dendrogram_name(self,field,axis):
    	self.get_fits_name(field,axis)
        self.dendrogram_name = "%s/%s_dendrogram.fits" %(self.fits_path, field)
        return self.dendrogram_name

    def get_fits_array(self,field,axis):
        self.get_fits_name(field,axis)
        if field == 'density':
            fits_fptr = pyfits.open(self.density_fits_name)  
            data = fits_fptr[0].data
            return data
            #return pyfits.open(self.density_fits_name)[0].data
        else:
            fits_fptr = pyfits.open(self.fits_name)  
            data = fits_fptr[0].data
            #data = np.array(fits_fptr[0].data).astype('double')
            return data
            #return pyfits.open(self.fits_name)[0].data
    
    def get_plot_location(self):
        # Create intermediate directories to store plots if they do not exist
        if not os.path.exists(self.plot_location):
            os.system("mkdir -p %s" %(self.plot_location))
        return self.plot_location

    def get_plot_location_P_frac(self):
        # Create intermediate directories to store plots if they do not exist
        if not os.path.exists(self.plot_location_P_frac):
            os.system("mkdir -p %s" %(self.plot_location_P_frac))
        return self.plot_location_P_frac        

    def get_plot_location_field_strength(self):
        if not os.path.exists(self.plot_location_field_strength):
            os.system("mkdir -p %s" %(self.plot_location_field_strength))
        return self.plot_location_field_strength

    def get_plot_location_viewing_angle(self):
        if not os.path.exists(self.plot_location_viewing_angle):
            os.system("mkdir -p %s" %(self.plot_location_viewing_angle))
        return self.plot_location_viewing_angle  
        
    def get_plot_location_time_dependence(self):
        if not os.path.exists(self.plot_location_time_dependence):
            os.system("mkdir -p %s" %(self.plot_location_time_dependence))
        return self.plot_location_time_dependence

    def get_plot_location_threshold_density(self):
        if not os.path.exists(self.plot_location_threshold_density):
            os.system("mkdir -p %s" %(self.plot_location_threshold_density))
        return self.plot_location_threshold_density       

    def get_plot_location_leaves(self):
        # Create intermediate directories to store plots if they do not exist
        if not os.path.exists(self.plot_location_leaves):
            os.system("mkdir -p %s" %(self.plot_location_leaves))
        return self.plot_location_leaves                       

    def get_plot_location_P_theta_vs_B_theta(self):
        # Create intermediate directories to store plots if they do not exist
        if not os.path.exists(self.plot_location_P_theta_vs_B_theta):
            os.system("mkdir -p %s" %(self.plot_location_P_theta_vs_B_theta))
        return self.plot_location_P_theta_vs_B_theta

    def get_plot_location_time_dependence(self):
        # Create intermediate directories to store plots if they do not exist
        if not os.path.exists(self.plot_location_time_dependence):
            os.system("mkdir -p %s" %(self.plot_location_time_dependence))
        return self.plot_location_time_dependence

    def get_plot_location_leaf_contours(self):
        # Create intermediate directories to store plots if they do not exist
        if not os.path.exists(self.plot_location_leaf_contours):
            os.system("mkdir -p %s" %(self.plot_location_leaf_contours))
        return self.plot_location_leaf_contours     

    def get_plot_location_mass_vs_smoothness(self):
        # Create intermediate directories to store plots if they do not exist
        if not os.path.exists(self.plot_location_mass_vs_smoothness):
            os.system("mkdir -p %s" %(self.plot_location_mass_vs_smoothness))
        return self.plot_location_mass_vs_smoothness  

    def get_plot_location_QU2EB(self):
        # Create intermediate directories to store plots if they do not exist
        if not os.path.exists(self.plot_location_QU2EB):
            os.system("mkdir -p %s" %(self.plot_location_QU2EB))
        return self.plot_location_QU2EB 
                 

    def get_data_location(self):
        # Create intermediate directories to store pickled data if they do not exist
        if not os.path.exists(self.data_location):
            os.system("mkdir -p %s" %(self.data_location))
        return self.data_location            

    def get_dendrogram(self, field, axis):
        self.get_dendrogram_name(field,axis)
        return Dendrogram.load_from(self.dendrogram_name)

    def make_frb(self,field, axis, res):
        self.get_fits_name(field,axis)
        print 'in make frb',self.fits_name 
        if field == 'density' and glob.glob(self.density_fits_name) != []:
            print "Existing file at the following directory.  I won't overwrite."
            print self.density_fits_name
            return -1
        if glob.glob(self.fits_name) != []:
            print "Existing file at the following directory.  I won't overwrite."
            print self.fits_name
            return -1

        ds = yt.load(self.enzo_dataset)

        stokes = re.compile('[QUI][xyz]')    # Regex matching stokes: 'Qx', 'Uz', etc.
        if re.match(stokes, field):
            # Use this for stokes
            #print "\n\nfield_parameters={'n0':%d, 'p':%d})\n"%(self.n0,self.p)
            #proj = ds.proj(field, axis, field_parameters={'n0':self.n0, 'p':self.p}) 
            proj = ds.proj(field, axis)                    
        if field == 'Bx' or field == 'By' or field == 'Bz':
            # Use this only for B field projections
            print '\n\n\nmaking', field, 'projection\n\n\n'
            proj = ds.proj(field, axis, weight_field='density')  
        else:
            proj = ds.proj(field,axis)

        # print "Writing yt projection"
        # print proj.to_pw().save()
        print '\n****** Making FRB ******\n'
        frb = proj.to_frb(1,res)
        hdu = pyfits.PrimaryHDU(frb[field])
        hdulist = pyfits.HDUList([hdu])

        # Create intermediate directories in fits_location if they do not exist
        if not os.path.exists(self.fits_path):
            os.system("mkdir -p %s" %(self.fits_path))

        if field == 'density':
            hdulist.writeto(self.density_fits_name)
            print "writing fits\n", self.density_fits_name
        else:
            hdulist.writeto(self.fits_name)
            print "writing fits\n", self.fits_name

    def make_projection(self, field, axis):
        ds = yt.load(self.enzo_dataset)
        return yt.ProjectionPlot(ds,axis,field)


    def make_dendogram(self, field, axis, res):
        """Make a dendrogram of a field on a specified axis
        and dump the dendrogram object to disk.
        """
        self.get_fits_name(field,axis)
        self.get_dendrogram_name(field,axis)
        if glob.glob(self.dendrogram_name) != []:
            print "Existing file at the following directory.  I won't overwrite."
            print self.dendrogram_name
            return -1
        print 'making dendrogram'
        data = self.get_fits_array(field, axis)
        d = Dendrogram.compute(data, min_value=2.0, min_delta=1.0, min_npix=5, verbose=True)
        d.save_to(self.dendrogram_name)
        return d

