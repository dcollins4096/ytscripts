
execfile('go')

#-----------------
#Gaussian Function
#----------------
def gauss(x, off,a, b):
    return off * np.exp(-(((x-a)**2)/(2*(b**2)))) 
def gauss_fit(centers,this_spectra):
    norm = this_spectra.sum()
    vbar= (this_spectra*centers).sum()/norm
    sigma2 = (this_spectra*(centers-vbar)**2).sum()/norm
    sigma = (sigma2)**0.5
    halfmax=this_spectra.max()*0.5
    #plt1.plot([vbar-0.5*sigma,vbar+0.5*sigma],[halfmax,halfmax])
    #ydata = gauss(centers, vbar, sigma)
    popt, pcov = curve_fit(gauss, centers, this_spectra, p0=[this_spectra.max(),vbar,sigma] ) #,  bounds=(0, [(1.1*np.max(this_spectra)), 0.15, np.max(this_spectra)]), method='trf')
    return {'vbar':vbar,'norm':norm,'sigma':sigma,'fit_norm':popt[0],'fit_center':popt[1],'fit_width':popt[2]}

from astrodendro import Dendrogram as dg
import pyfits
import fPickle
import clfind2d
from yt.analysis_modules.ppv_cube.api import PPVCube
import yt.units as u
import scipy as sp
import scipy.optimize as opt
from scipy.optimize import curve_fit
import clump_properties
reload(clump_properties)
import mw_stuff
import clump_subset
reload(mw_stuff)
class access_thing():
    def __init__(self, sim, frame, machine='nazare',frb_resolution=None,fits_basedir=None):
        """Wraps up the names and locations of things for ease of use and platform interchangability.
        *sim* can be: high_512, mid_512, low_512, test_uniform, test_ot
        We need to get out the location of the enzo data,
        and the location of the related fits projection.
        """
        self.sim=sim
        self.frame=frame
        restart_sim_template = "%s/RS%04d/restart%04d"
        dataset_sim_template = "%s/DD%04d/data%04d"

        """start spaghetti of names and locations"""

        self.Paper08_sims = ['high_512','mid_512','low_512','high_256','mid_256','low_256']
#       self.fits_location_dict = {'low_512':'Low_512/fit_b205v','mid_512':'Mid_512/fit_b25v', 'high_512':'High_512/fit_b025v',
#                                  'test_uniform':'TEST_uniform/fits/fits','test_ot':'/TestOT/fits/fits',
#                                  'test_sphere':'test_sphere/fits'}

        self.fits_location_dict={'high_512':'/B02/512','high_256':'/B02/256',
                                 'low_256':'/B20/256', 'low_512':'/B20/512'}
        self.sim_dir = {'high_512':'/B02/512','mid_512':'/B2/512','low_512':'/B20/512',
                        'high_256':'/B02/256','mid_256':'/B2/256','low_256':'/B20/256',
                        'test_uniform':'TEST_uniform/run',
                        'test_ot':'TestOT/run',
                        'test_sphere':'test_sphere/run_link'}[sim]
        if machine == 'nazare':
            if sim in self.Paper08_sims:
                self.enzo_basedir = '/scratch1/dcollins/Paper08'
            if sim in ['test_uniform','test_ot', 'test_sphere']:
                self.enzo_basedir = '/scratch1/dcollins/Paper41_rht'
        if fits_basedir is not None:
            self.fits_basedir = '/scratch1/dcollins/Paper14b/'
        else:
            self.fits_basedir = fits_basedir
            


        if sim in self.Paper08_sims:
            self.enzo_dataset_template = restart_sim_template
        else:
            self.enzo_dataset_template = dataset_sim_template
            
        """ end spaghetti """


        self.enzo_location = "%s/%s"%(self.enzo_basedir,self.sim_dir)
        self.fits_location = self.fits_location_dict[sim]
        
        self.enzo_dataset = self.enzo_dataset_template%(self.enzo_location, frame,frame)
        self.resx=256
        self.resy=256
        if frb_resolution is not None:
            self.resx = frb_resolution
            self.resy = frb_resolution

        self.ds = None
        self.cut_region = None
        self.leaf_clumps = {}
        self.clump_property_list = {}
        self.master_clumps = {}
        self.cut_regions = {}
        self.clump_sets = {}

    def get_fits_name(self,field,axis):
        self.fits_dir = "%s/%s/%s/RS%04d"%(self.fits_basedir,self.fits_location,axis,self.frame)
        subsubdir=''
        for subdir in self.fits_dir.split('/'):
            if subdir is not '':
                subsubdir +='/'+ subdir
                glb =  glob.glob(subsubdir)
                if glb == []:
                    os.mkdir(subsubdir)
        self.fits_name = "%s/%s.fits"%(self.fits_dir,field)
        
    def get_fits_array(self,field,axis):
        self.get_fits_name(field,axis)
        if glob.glob(self.fits_name) == []:
            self.make_frb(field,axis)
        return pyfits.open(self.fits_name)[0].data

    
    def make_frb(self,field, axis):
        self.get_fits_name(field,axis)
        if glob.glob(self.fits_name) != []:
            print "Existing file at the following directory.  I won't overwrite."
            print self.fits_name
            return -1
        ds = yt.load(self.enzo_dataset)
        proj = ds.proj(field,axis,center=[0.5]*3)
        print "Writing yt projection"
        print proj.to_pw().save()
        frb = proj.to_frb(1,[self.resx,self.resy])
        hdu = pyfits.PrimaryHDU(frb[field])
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(self.fits_name)
        print "writing fits\n", self.fits_name

    def get_dendrogram(self,axis,prefix='dendrogram_1',fpickle_args=None):
        """Unsure how to save these to disk"""
        self.get_fits_name('density',axis)
        dendro_filename = self.fits_dir + "/%s.fits"%prefix
        if len(glob.glob(dendro_filename)) > 0:
            pyfits.open(dendro_filename)[0].data
            #return fPickle.load(dendro_filename)
        else:
            density = self.get_fits_array('density',axis)
            if fpickle_args is None:
                fpickle_args={'min_value':2.0, 'min_delta':1.0, 'min_npix':10.0, 'verbose':True}
            d = dg.compute(density, **fpickle_args)
            
            #fPickle.dump(d,dendro_filename)
            return d
    def get_cl2d(self,axis,prefix='cl2d_1',clfind_args=None):
        density = self.get_fits_array('density',axis)
        name_root = '%s/%s'%(self.fits_dir,prefix)
        outname = "%s_Mask.fits"%name_root
        if clfind_args is None:
            clfind_args = {'levels':range(1,12), 'log':True,'nPixMin':3}
        if len(glob.glob(outname)) == 0:
            clfind2d.clfind2d(self.fits_name,name_root,**clfind_args)
        return pyfits.open(outname)[0].data

    def image_cl2d(self,axis,clump_prefix='CLUMP',ppv_prefix=None,out_prefix='p14b_image1', plot_number=True):
        if ppv_prefix:
            cube_read,cube_ds,vel = self.get_ppv(axis,ppv_prefix)
            density_ppv = np.sum(cube_read,axis=2)
        density_fits = self.get_fits_array('density',axis)
        clumps = self.get_cl2d(axis,clump_prefix)
        nx,ny = clumps.shape
        fig = plt.figure()
        plt1 = fig.add_subplot(111)
        plt1.imshow(np.log10(density_fits),origin='lower',interpolation='nearest',cmap='gray')
        clumps_contour_wtf = 1.0*copy.copy(clumps)/clumps.max()
        clumps_contour_wtf[clumps >0 ] = 1
        plt1.contour(clumps_contour_wtf)
        if plot_number:
            y,x = np.mgrid[0:ny,0:nx]
            for n in np.unique(clumps)[1:]:
                n_x=np.mean(x[clumps == n])
                n_y=np.mean(y[clumps == n])
                plt1.text(n_x,n_y,"%d"%n, color='g')
                print "centroid", n_x,n_y

        outname = out_prefix+"_fits_clumps.png"
        fig.savefig(outname)

        plt.close(fig)
        if ppv_prefix:
            plt.clf()
            plt.imshow(np.log10(density_ppv), origin='lower',interpolation='nearest',cmap='gray')
            outname_image = out_prefix+"_ppv.png"
            plt.savefig(outname_image)

        plt.clf()
        plt.imshow(clumps, origin='lower',interpolation='nearest')
        outname_clumps = out_prefix+"_clumps.png"
        plt.savefig(outname_clumps)
        if ppv_prefix:
            print "wrote: PPV", outname_image
        print "wrote: FITS FRB", outname
        print "wrote: clumps", outname_clumps

    def get_ppv(self,axis,prefix='ppv_1',ppv_args={}):
        self.get_fits_name('density',axis)
        ds = yt.load(self.enzo_dataset)
        ppv_name_vel = '%s/%s_vel.fits'%(self.fits_dir,prefix)
        ppv_name_yt = '%s/%s_yt.fits'%(self.fits_dir,prefix)
        normal = axis
        if len(glob.glob(ppv_name_yt)) == 0:
            field = ppv_args.get('field','density')
            velocity_bounds = ppv_args.get('velocity_bounds',(-35,35,14,'code_velocity'))
            dims = ppv_args.get('dims',256)
            cube = PPVCube(ds, normal, "density", velocity_bounds,dims=dims, method='sum')
            cube.write_fits(ppv_name_yt)
            hdu = pyfits.PrimaryHDU(cube.vbins)
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto(ppv_name_vel,clobber=True)
            print "Wrote", ppv_name_yt
            print "Wrote", ppv_name_vel
        print "Reading", ppv_name_yt
        #self.ppv_name = ppv_name
        self.ppv_name_yt = ppv_name_yt
        ds_cube = yt.load(ppv_name_yt)
        cube_only = np.rot90(pyfits.open(ppv_name_yt)[0].data.swapaxes(0,2),3)
        vel =  pyfits.open(ppv_name_vel)[0].data
        return cube_only, ds_cube, vel

    def masses(self,axis,prefix,background_column=0):
        mask = self.get_cl2d(axis,prefix)
        density = self.get_fits_array('density',axis)
        masses = np.zeros(mask.max()+1)
        for n in range(1,mask.max()+1):
            this_density_list = density[mask == n] - background_column
            masses[n]= this_density_list.sum()/(self.resx*self.resy) 
        return nar(masses)

    def spectra(self,axis,clump_prefix,ppv_prefix):
        mask = self.get_cl2d(axis,clump_prefix)
        cube_read,cube_ds,vel = self.get_ppv(axis,ppv_prefix)
        nx,ny,nz = cube_read.shape
        spectra = np.zeros([mask.max()+1,nz])
        flat_cube = np.sum(cube_read,axis=2)
        for clump in np.unique(mask)[1:]: #[1:]: #kludge was 1:
            this_mask = mask == clump
            for n in range(nz):
                this_2d = cube_read[:,:,n][this_mask]
                spectra[clump,n] = this_2d.sum()
        return vel,spectra

    def spectra_image_pair(self, axis,clump_prefix,ppv_prefix,image_zoom=True):
        density = self.get_fits_array('density',axis)
        mask = self.get_cl2d(axis,clump_prefix)
        vel, spectra = self.spectra(axis,clump_prefix,ppv_prefix)
        centers = 0.5*(vel[1:]+vel[:-1])
        cube_read,cube_ds,vel = self.get_ppv(axis,ppv_prefix)
        nx,ny,nz = cube_read.shape
        y,x = np.mgrid[0:ny,0:nx]
        for clump in np.unique(mask)[1:]: #[1:]: #kludge was 1:
            this_spectra = spectra[clump]
            if np.isnan(this_spectra).sum() > 0:
                print "Skipping clump", clump, "contains nans"
                continue
            fig = plt.figure()
            plt1 = fig.add_subplot(121)
            plt1.plot(centers,this_spectra)
            fit = gauss_fit(centers,this_spectra)
            plt1.plot(centers,gauss(centers,fit['fit_norm'],fit['fit_center'],fit['fit_width']),c='r')
            
            plt2 = fig.add_subplot(122)
            x_slice = slice(x[mask==clump].min(),x[mask==clump].max()+1)
            y_slice = slice(y[mask==clump].min(),y[mask==clump].max()+1)
            if image_zoom:
                #Gaah, I'm having a stroke, why is this the order of slices?
                plt2.imshow( np.log10(density[y_slice,x_slice]),origin='lower',interpolation='nearest',cmap='gray')
            else:
                density_to_plot = copy.copy(density)
                density_to_plot[mask != clump] = density_to_plot[mask == clump].min()
                plt2.imshow( np.log10(density_to_plot),origin='lower',interpolation='nearest',cmap='gray')
            outname = 'p14b_2bpanel_clump_%04d.png'%clump
            fig.savefig(outname)
            print outname
            plt.close(fig)
        return vel,spectra

    def alpha(self,axis,clump_prefix,ppv_prefix, background_column=0):
        vel,spectra = self.spectra(axis,clump_prefix,ppv_prefix)
        centers = 0.5*(vel[1:]+vel[:-1])
        mask = self.get_cl2d(axis,clump_prefix)
        masses =self.masses(axis,clump_prefix,background_column)
        ds = yt.load(self.enzo_dataset)
        G = ds['GravitationalConstant']/(np.pi*4.)
        alpha = np.zeros(mask.max()+1)
        for clump in np.unique(mask)[1:]: #[1:]: #kludge was 1:
            this_size = ((mask==clump).sum()*1./(self.resx*self.resy))**0.5
            this_spectra = spectra[clump]
            if np.isnan(this_spectra).sum() > 0:
                print "Skipping clump", clump, "contains nans"
                alpha[clump]=0
                continue
            fit = gauss_fit(centers,this_spectra)
            sigma2 = fit['fit_width']**2
            alpha[clump]= 5.*sigma2*this_size/(masses[clump]*G)
        return alpha

    def alpha_mass(self,axis,clump_prefix,ppv_prefix, background_column=0, with_numbers=False):
        alpha = self.alpha( axis,clump_prefix,ppv_prefix,background_column)
        masses = self.masses(axis,clump_prefix,background_column)
        plt.clf()
        plt.xscale('log')
        plt.yscale('log')
        #ef('p14_data.py')
        import p14_data
        M1,A1 = p14_data.get_real_data()
        SolarMassesPerCodeMass = 5900
        plt.scatter( M1, A1, c='k', marker=7, label = 'Real Data')
        non_zero = alpha > 0
        if with_numbers:
            for n in range(len(masses)):
                if alpha[n] > 0:
                    plt.text(masses[n]*SolarMassesPerCodeMass, alpha[n], '%d'%n)
        else:
            plt.scatter(masses[non_zero]*SolarMassesPerCodeMass,alpha[non_zero])
        x_min = min([M1.min(), (masses[non_zero]*SolarMassesPerCodeMass).min()])
        x_max = max([M1.max(), (masses[non_zero]*SolarMassesPerCodeMass).max()])
        plt.xlim(x_min*0.9,x_max*1.1)
        y_min = min([A1.min(), alpha[non_zero].min()])
        y_max = max([A1.max(), alpha[non_zero].max()])
        plt.ylim(y_min*0.9,y_max*1.1)
        plt.xlabel('M')
        plt.ylabel('alpha')
        plt.savefig('p14b_malpha.png')

    def get_cut_region(self,axis,clump_prefix,clump_index):
        self.twod_clump = np.transpose(self.get_cl2d(axis,clump_prefix))
        flat_clump = self.twod_clump.flatten()
        clump_stuff_tuple = (flat_clump,'xyz'.index(axis),self.twod_clump)

        nx, ny = self.twod_clump.shape
        ind_map = ((na.mgrid[0:nx, 0:ny] + 0.5)/nx) 

        these_x = ind_map[0][self.twod_clump == clump_index].flatten()
        these_y = ind_map[1][self.twod_clump == clump_index].flatten()
        Left = na.zeros(3)
        Right = na.ones(3)
        Left[ x_dict[0] ]  = na.floor(these_x.min()*nx)/nx
        Left[ y_dict[0] ]  = na.floor(these_y.min()*nx)/nx
        Right[ x_dict[0] ] = na.ceil(these_x.max()*ny)/ny
        Right[ y_dict[0] ] = na.ceil(these_y.max()*ny)/ny 
        Center = 0.5*(Left+Right)                         
        if self.ds is None:
            self.ds = yt.load(self.enzo_dataset)
        region = self.ds.region(Center,Left,Right)
        region.set_field_parameter('clump_mask_stuff',clump_stuff_tuple)
        self.cut_region = region.cut_region(['obj["clump_mask"].astype("int") == %d'%clump_index])
        self.cut_region.set_field_parameter('clump_mask_stuff',clump_stuff_tuple)


    def get_clumps(self,axis,clump_prefix,clump_index, clump_parameters={}):
        self.get_cut_region(axis,clump_prefix,clump_index)
        master_clump = Clump(self.cut_region,"density")
        master_clump.add_validator("min_cells", 20)
        clump_parameters['c_min'] = clump_parameters.get('c_min',10)
        clump_parameters['c_max'] = clump_parameters.get('c_max',self.cut_region['density'].max())
        clump_parameters['step'] = clump_parameters.get('step',2)
        find_clumps(master_clump, clump_parameters['c_min'], clump_parameters['c_max'], clump_parameters['step'])
        self.leaf_clumps[clump_index] = get_lowest_clumps(master_clump) #if both min_cells and grav_bound are used, this is empty.
        self.clump_index = clump_index
        self.working_clump_prefix = clump_prefix
        self.master_clumps[clump_index] = master_clump
        self.cut_regions[clump_index] = self.cut_region
        return self.leaf_clumps[clump_index]
    
    def to_mw(self):
        if self.leaf_clumps.has_key[self.clump_index] is None:
            for n_sub,clump in enumerate(self.leaf_clumps[self.clump_index]):
                #directory = "/lustre/medusa/collins/Paper12/MWgrav/simplegrav"
                directory = "%s/clumps_grav_%s"%(self.fits_dir,self.working_clump_prefix)
                glb =  glob.glob(directory)
                if glb == []:
                    os.mkdir(directory)
                name_base = "%s/mw_clump_cl%04d_sub%04d"%(directory,self.clump_index,n_sub)
                mw_stuff.run_mw_grav(name_base,  clump)
        else:
            print "make clumps first"



    def plot_clumps(self, output_prefix = 'p14b_thing'):
        if self.leaf_clumps[self.clump_index] is None:
            print "run clumps please"
        else:
            if self.ds is None:
                self.ds = yt.load(self.enzo_dataset)
            for ax in [0,1,2]:
                proj = yt.ProjectionPlot(self.ds,ax,'density')
                proj.set_cmap('density','gray')
                proj.annotate_clumps(self.leaf_clumps[self.clump_index])
                print proj.save(output_prefix)

    def load_clumps(self,axis,clump_prefix,clump_index):
        self.get_fits_name('density','x')
        ds = yt.load(self.enzo_dataset)
        self.working_clump_prefix = clump_prefix
        self.clump_index = clump_index
        clump_dir = "%s/clumps_%s"%(self.fits_dir,self.working_clump_prefix)
        clump_glob = "%s/clump%04d_sub????.pickle"%(clump_dir,self.clump_index)
        self.leaf_clumps[clump_index] = []
        self.ds_list = []
        fname_list = sorted(glob.glob(clump_glob))
        for  fname in fname_list[0:3]:
            print "kludge"
            this_ds, this_clump = fPickle.load(fname)
            self.leaf_clumps[clump_index].append(this_clump)
            #self.leaf_clumps.append(Clump(this_clump,'density'))
            self.ds_list.append(this_ds)
            print fname

    def compute_clump_properties(self):
        if self.leaf_clumps.has_key(self.clump_index):
            self.clump_dir = "%s/clumps_%s"%(self.fits_dir,self.working_clump_prefix)
            clump_stuff_fname = "%s/clump_stuff.pickle"%self.clump_dir
            self.clump_property_list[self.clump_index]=[]
            if len(glob.glob(clump_stuff_fname) ) > 0 and False:
                clump_prop_list_temp = fPickle.load(clump_stuff_fname)
                for clump_ind_tmp in clump_prop_list_temp:
                    self.clump_property_list[clump_ind_tmp] = clump_prop_list_temp[clump_ind_tmp]
            if len(glob.glob(self.clump_dir)) == 0:
                os.mkdir(self.clump_dir)

            for n,L in enumerate(self.leaf_clumps[self.clump_index]):
                self.clump_property_list[self.clump_index].append(clump_properties.clump_stuff(L,self.ds, 1,self.clump_index,n))


            fPickle.dump(self.clump_property_list,clump_stuff_fname)
            self.clump_sets[self.clump_index] = clump_subset.subset(self.clump_property_list[self.clump_index])
            print "dumped", clump_stuff_fname



    def save_clumps(self):
        if self.leaf_clumps.has_key(self.clump_index) is None:
            print "run clumps please"
        else:
            self.clump_dir = "%s/clumps_%s"%(self.fits_dir,self.working_clump_prefix)
            glb =  glob.glob(self.clump_dir)
            if glb == []:
                os.mkdir(self.clump_dir)
            for n,L in enumerate(self.leaf_clumps):
                fname = "%s/clump%04d_sub%04d.pickle"%(self.clump_dir,self.clump_index,n)
                print fname
                fPickle.dump(L.data,fname)
    def pdf(self,axis):
        density = self.get_fits_array('density',axis)
        plt.clf()
        plt.hist(np.log10(density.flatten()), histtype='step',bins=100)
        plt.yscale('log')
        outname = 'column_density.pdf'
        plt.savefig(outname)
        print outname



