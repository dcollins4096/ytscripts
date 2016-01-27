import pyfits
import scatter_fit
reload(scatter_fit)

def log_fun(field):
    return np.log10(field)

if 'use_high_res' not in dir():
    use_high_res=True
error_list = []

if 'fig_list' in dir():
    for fig in fig_list:
        print 'delfig'
        plt.close(fig)
    del fig_list
fig_list=[]
def newfig():
    print "NEWFIG"
    fig_list.append(plt.figure())
    return fig_list[-1]
def point_selector(x,y,x0,x1,y0,y1,dx,dy):
    #F(x,y) = (y1-y0)*x + (x0-x1)*y +(x1*y0-x0*y1)
    #F==0, (x,y) is on the line.  <0 above, >0 below
    #All above or below, then there's no intersection. 
    #Sign change means that two points are on opposite sides of the line.
    c1 = y1-y0
    c2 = x0-x1
    c3 = x1*y0-x0*y1
    a   = c1*(x      )+c2*(y      )+c3
    amm = c1*(x-dx/2.)+c2*(y-dy/2.)+c3
    amp = c1*(x-dx/2.)+c2*(y+dy/2.)+c3
    app = c1*(x+dx/2.)+c2*(y+dy/2.)+c3
    apm = c1*(x+dx/2.)+c2*(y-dy/2.)+c3
    keep = amm*amp <= 0
    keep = np.logical_or( keep, amp*app <= 0 )
    keep = np.logical_or( keep, app*apm <= 0 )
    keep = np.logical_or( keep, apm*amm <= 0 ) 
    #pdb.set_trace()
    return keep

class filament():
    def __init__(self,nfilament, x_fil=None,y_fil=None, x_coords=None,y_coords=None,parent=None):
        """filament objects take a set of points in x and y, finds the centroid.  Really just a container."""
        self.x_fil=x_fil
        self.y_fil=y_fil
        if nfilament in [48]:
            #some filaments we want to curate.  I don't have a particularly good way to do this right now.
            keep = x_fil < 15./256
            self.x_fil = self.x_fil[keep]
            self.y_fil = self.y_fil[keep]
            
        """Check for periodic jumps."""
        if (x_fil < parent.dx[0]).any() or (y_fil < parent.dx[1]).any():
            print "Filament", nfil, "probably has periodic wrap, its against the edge of the domain."
        self.x_centroid = self.x_fil.sum()/self.x_fil.size
        self.y_centroid = self.y_fil.sum()/self.y_fil.size
        self.n_points = self.x_fil.size
        self.parent=parent


class perpandicular():
    def __init__(self,fil, ind, width=0.3, n_points_to_fit = 4):
        #fit = scatter_fit.scatter_fit(None,x_fil,y_fil, plot_points=False)#x_fil,y_fil)
        #slope = fit['fit'][0]
        #offset = fit['fit'][1]
        #ind = x_fil.size/2 
        spine_point_x = fil.x_fil[ind]; spine_point_y = fil.y_fil[ind]
        """the transverse line"""
        r = (spine_point_x-fil.x_fil)**2 + (spine_point_y - fil.y_fil)**2
        r_args = np.argsort(r)
        x_closest = fil.x_fil[r_args][0:n_points_to_fit]
        y_closest = fil.y_fil[r_args][0:n_points_to_fit]
        #poly fit has a problem with vertical lines.
        self.vertical_problem = np.abs((x_closest-x_closest.mean())).max() < 0.5*fil.parent.dx[0]
        self.horizontal_problem = np.abs((y_closest-y_closest.mean())).max() < 0.5*fil.parent.dx[1]
        if self.vertical_problem:
            Deltax = width
            Deltay = 0
            #y0 = max([spine_point_y - 0.5*Deltay,0.5]); y1 = min([spine_point_y + 0.5*Deltay, y.max()])
            self.y0=self.y1=spine_point_y
            self.x0 = max([spine_point_x - 0.5*Deltax,fil.parent.dx[0]]); self.x1 = min([spine_point_x + 0.5*Deltax, fil.parent.x.max()])
        else:
            #fit = scatter_fit.scatter_fit(None,x_closest,y_closest, plot_points=False)#x_fil,y_fil)
            to_plot = None; #low_res_ax
            fit = scatter_fit.scatter_fit(to_plot,x_closest,y_closest, plot_points=False)#x_fil,y_fil)
            slope = fit['fit'][0]
            offset = fit['fit'][1]
            perp_slope = -1./slope
            Theta = np.arctan(perp_slope)
            if slope < 1e-5:
                self.horizontal_problem = True
            #Spatial extent in line in X and Y
            Deltax = width * np.cos(Theta) 
            Deltay = width * np.sin(Theta) 
            self.x0 = max([spine_point_x - 0.5*Deltax,fil.parent.dx[0]]); self.x1 = min([spine_point_x + 0.5*Deltax, fil.parent.x.max()])
            self.y0 = max([spine_point_y - 0.5*Deltay,fil.parent.dx[1]]); self.y1 = min([spine_point_y + 0.5*Deltay, fil.parent.y.max()])
        #for the bounding box, we want the ordered set.  The bounding box will be used to extract subsets of data
        self.yLeft  = min([self.y0,self.y1]);  
        self.yRight = max([self.y0,self.y1]);  
        self.xLeft  = min([self.x0,self.x1]);  
        self.xRight = max([self.x0,self.x1]);  
        self.spine_point_x = spine_point_x
        self.spine_point_y = spine_point_y

class filament_tool():
    """Main image container object."""
    def __init__(self, source_file,frame=70, resolution=[256,256],box_size_pc = 4.6,
                extraction_width_pc=0.3):
        self.filament_dict={}
        self.setname = setname
        self.filename_prefix = filename_prefix
        self.resolution = np.array(resolution)
        self.data_nobuf = pyfits.open(source_file)[0].data
        self.data = np.zeros(resolution)
        self.frame=frame
        resolution_offset = [ (R-N)/2 for R,N in zip(resolution, self.data_nobuf.shape)]
        print "RES OFF", resolution_offset
        if resolution_offset[0]*resolution_offset[1] > 0:
            self.data[resolution_offset[0]:-resolution_offset[0],resolution_offset[1]:-resolution_offset[1]] = self.data_nobuf
        else:
            self.data = self.data_nobuf
        self.dx = 1./self.resolution
        self.Right = np.ones(2)
        self.y,self.x = np.mgrid[0.5*self.dx[0]:self.Right[0]-0.5*self.dx[0]:self.resolution[0]*1j,
                                 0.5*self.dx[1]:self.Right[0]-0.5*self.dx[1]:self.resolution[1]*1j]
        self.box_size_pc = box_size_pc
        self.extraction_width_pc = extraction_width_pc
        self.profile_aggregator = None
        self.make_fig()
    def make_fig(self):
        self.image_fig = newfig()
        self.image_ax = self.image_fig.add_subplot(111)
        self.profile_fig = newfig()
        self.profile_ax = self.profile_fig.add_subplot(111)
    def make_image(self,cmap='gray', log=False):
        if log:
            to_plot = np.log10(self.data)
        else:
            to_plot = self.data
        d1=self.image_ax.imshow(to_plot, interpolation='nearest',origin='lower',cmap=cmap,extent=[0,1,0,1]) #extents are left,right,bottom,top
        self.image_fig.colorbar(d1)
        self.profile_aggregator = None
    def image_save(self,outname):
        self.image_ax.set_xlim(0,1)
        self.image_ax.set_ylim(0,1)

        self.image_fig.savefig(outname)
        print outname
    def profile_save(self,outname, take_log=False):
        #outname = 'filament_image_%s_n%04d_r%04d.png'%(self.setname, self.frame, self.resolution)
        if take_log:
            self.profile_fig.set_yscale('log')
        self.profile_ax.set_ylabel('density')
        self.profile_ax.set_xlabel('x[pc]')

        self.profile_fig.savefig(outname)
        print outname
    def plot_filament(self,this_filament, c='b', plot_number=None):
        self.image_ax.scatter(this_filament.x_fil, this_filament.y_fil, marker='o', c=c, linewidths=0,s=1)
        if plot_number is not None:
            self.image_ax.text(this_filament.x_centroid,this_filament.y_centroid, "%d"%plot_number, fontsize=5, color=c)

    def get_filament_points(self,nfilament):
        """subclassing should start here."""
        filament_mask = self.data == nfil
        x_fil = self.x[ filament_mask ]
        y_fil = self.y[ filament_mask ]
        self.filament_dict[nfilament] = filament(nfilament,x_fil,y_fil, parent=self)

        return self.filament_dict[nfilament]

    def extract_profile(self, perp, shift_peak=True):
        #self.image_ax.plot([perp.xLeft,perp.xRight], [perp.yLeft,perp.yRight], c='g')
#        self.image_ax.plot([perp.x0,perp.x1], [perp.y0,perp.y1], c='y')
        self.image_ax.scatter([perp.spine_point_x],[perp.spine_point_y], c='g')
        slice_x = slice(np.where(self.x[0,:]<=perp.xLeft)[0][-1], np.where(self.x[0,:]>=perp.xRight)[0][0]+1)
        slice_y = slice(np.where(self.y[:,0]<=perp.yLeft)[0][-1], np.where(self.y[:,0]>=perp.yRight)[0][0]+1)
        x_sub=self.x[slice_y,slice_x]
        y_sub=self.y[slice_y,slice_x]
        keep = point_selector(x_sub,y_sub,perp.x0,perp.x1,perp.y0,perp.y1,self.dx[0],self.dx[1])
        x_a = x_sub[keep]
        y_a = y_sub[keep]
        density = self.data[slice_y,slice_x][keep] # fullsetnar([fullset[ix,iy] for ix,iy in zip(x_i,y_i)])
        new_center_x = perp.spine_point_x; new_center_y = perp.spine_point_y
        if shift_peak:
            """ Shift the curve to match the peak.  To avoid multiple maxima, restrict to a window of 0.1 pc """
            peak_hunt_subset = (x_a - new_center_x)**2+(y_a-new_center_y)**2 < (0.05/4.6)**2 #0.05 pc in pixels
            max_coord=int( np.where(density[peak_hunt_subset]==density[peak_hunt_subset].max())[0].mean() )
            new_center_x =  x_a[peak_hunt_subset][max_coord]
            new_center_y =  y_a[peak_hunt_subset][max_coord]
        if (new_center_x-perp.spine_point_x)**2+(new_center_y-perp.spine_point_y)**2 > self.extraction_width_pc**2:
            error = "Error: Peak for filament %d index %d shift by more than 0.25 width"%(nfil, ind)
            error_list.append(error)
            print error
        self.image_ax.scatter(x_a,y_a, marker='o',c='r', linewidths=0, s=0.1)
        sign_of_line   = np.sign(x_a-new_center_x)
        sign_of_line_y = np.sign(y_a-new_center_y) #this funny sign juggle takes care of vertical points
        sign_of_line[sign_of_line==0] = sign_of_line_y[sign_of_line==0]
        coordinates = sign_of_line*np.sqrt((x_a-new_center_x)**2+(y_a-new_center_y)**2) #This centers the profile on the Filament.
        coordinates *= self.box_size_pc
        sort_coord = np.argsort(coordinates)
        self.profile_ax.plot(coordinates[sort_coord], density[sort_coord], marker='o',c=rmap(ind))
        if self.profile_aggregator is None:
            self.profile_aggregator ={ 'coord':[], 'density':[], 'all_coord':[] }
        self.profile_aggregator['coord'].append(coordinates)
        self.profile_aggregator['all_coord'] += coordinates.tolist()
        self.profile_aggregator['density'].append(density)

    def hist_profile(self):
        plt.clf()
        c=np.array(self.profile_aggregator['coord']).flatten()
        plt.hist(c,bins=100)
        plt.savefig('tmp.png')




"""these sets lost two zones to smoothing.  Restore."""
directory = '/Users/dcollins/RESEARCH2/Paper37_Philaments/2015-06-12-disperse/DATA'
filename_prefix = 'b02_512_0070_smoothed_0256_density'
frame = 70
setname = 'b02'

filament_map = filament_tool(source_file = '%s/%s.%s'%(directory,filename_prefix,'fits_c1.up.NDskl.fits'),
                             resolution = [256,256], frame=frame)

low_map = filament_tool('%s/%s_512_%04d_smoothed_%04d_density.fits'%(directory,setname,frame,256),
                        resolution = [256,256], frame=frame)
high_map = filament_tool(source_file= '%s/%s_512_%04d_%04d_projection_density.fits'%(directory,setname,frame,8192),
                         resolution = [8192]*2, frame=frame)


#filament_map.make_image(cmap='jet')
#filament_map.image_save('x2_filament_image_%s_n%04d_r%04d.png'%(filament_map.setname, filament_map.frame, filament_map.resolution[0]))
high_map.make_image(log=True)
#high_map.image_save('x2_high_res_no_filaments_%s_n%04d_r%04d.png'%(high_map.setname, high_map.frame, high_map.resolution[0]))
low_map.make_image(log=True)
#low_map.image_save('test.png')



filament_cmap=rainbow_map(filament_map.data.max()+1)
for nfil in np.unique(filament_map.data)[1:]:
    if nfil not in [42]:
        continue
    this_color = filament_cmap(nfil)
    this_one = filament_map.get_filament_points(nfil)
    low_map.plot_filament(this_one, c=this_color,plot_number=None ) #nfil)
    low_map.image_save('x2_low_f%04d.png'%nfil)
    high_map.plot_filament(this_one, c=this_color,plot_number=nfil)
    high_map.image_save('x2_high_f%04d.png'%nfil)

    """For actual profiles along the filament"""
    rmap = rainbow_map(this_one.n_points+1)
    for ind in  range(this_one.n_points):
        this_perp = perpandicular(this_one,ind,0.3/4.6,n_points_to_fit=-1)
        low_map.extract_profile(this_perp)
        #high_map.extract_profile(this_perp)
    low_map.image_save('low_f%04d.pdf'%nfil)
    low_map.profile_save('low_f%04d_profile_test.pdf'%nfil)
    #high_map.image_save('high_f%04d.pdf'%nfil)
    #high_map.profile_save('high_f%04d_profile.pdf'%nfil)

            

