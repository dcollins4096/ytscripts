#from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib import use; use('Agg') 
import time
import re
from astropy.io import fits
import matplotlib.pyplot as plt
import sys

sys.path.append("/home/huffenbe/projects/cmbtools")
#import cmbtools
import cmbtools_handler as cmbtools

start = time.time()

machine = 'nazare'
#machine = 'corey_local'

ef=execfile
ef('access_thing.py')
ef('fPickle.py')

res = 8192
#frame = 50
unweighted = False
theta_vs_r_whole_box = True
sims = ['High_512', 'Mid_512', 'Low_512']
axes = ['x','y','z']
params = [[19,0], [39,0], [1945,0], [1,1]]
frames = [10, 30, 50]
bins = range(0,91,10)

# Diameter of 0.1 pc (large) and 0.01 pc (small) disks for 8192 res
# or 1 and 0.1 for 512 res 
BOXSIZE_PC = 4.6                 # Simulation is 4.6 parsecs on one side
PC_TO_PIXEL = res/BOXSIZE_PC     # [pixels/parsec]
if res == 8192:
    r_sm_pc = 0.005
    r_lg_pc_max = 0.05            
elif res == 512:
    r_sm_pc = 0.05
    r_lg_pc_max = 0.5
else:
    print '**** res should be 8192 or 512 for 512^3 datasets ****'
    
# Convert disk sized in pc to integer sizes in pixels
r_sm = int(round(r_sm_pc * PC_TO_PIXEL))
r_lg_max = int(round(r_lg_pc_max * PC_TO_PIXEL))
    
# Make list of disk radii
if theta_vs_r_whole_box:    
    r_disk = np.array([r_sm_pc, 0.05, 0.1, 0.25, 0.5, 1., 2.])
else:
    num_lg_disks = 25     # Number of disks from small scale to large scale
    r_disk = np.linspace(r_sm_pc, r_lg_pc_max, num_lg_disks)

# change [0.05, 0.1, 0.25, 0.5, 1., 2.] to [005,010,025,050,100,200] for use in file names
r_disk_str = ['%03d'%(num*100) for num in r_disk]


def make_frbs():
    for frame in [frames[2]]:
        #for sim in sims:
        for sim in ['High_512']:
            beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]
        
            # param has the form [n0,p]
            #for param in params:
            for param in [[1,1]]:
        
                # Make fixed resolution buffers
                if 0:
                    at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                
                    for axis in axes:
                        field_horizontal = {'x':'By','y':'Bz','z':'Bx'}[axis]
                        field_vertical   = {'x':'Bz','y':'Bx','z':'By'}[axis]              	
                        at.make_frb('Q%s_n0-%04d_p-%d'%(axis,param[0],param[1]),axis,res)
                        at.make_frb('U%s_n0-%04d_p-%d'%(axis,param[0],param[1]),axis,res)
                        at.make_frb('N2%s_n0-%04d_p-%d'%(axis,param[0],param[1]),axis,res)
                        at.make_frb(field_horizontal, axis, res)
                        at.make_frb(field_vertical, axis, res)                    
                        at.make_frb('density', axis, res)
                        at.make_dendogram('density', axis, res)
                

def analyze_dense_regions(theta_vs_r_whole_box=True):                      
    """Calculate polarization angle and polarization fraction in in 
    large (0.1 pc) and small (0.01 pc) disks around each spot of high
    density.
    dg  : denrogram object
    res : resolution of simulation
    Q, U: stokes parameters
    N : column density 
    N2 : correction to column density because of inclination of B field in plane of sky
    B_hor: B field in horizontal direction
    B_vir: B field in vertical direction
    """
    pFrac_max = 0
    total_masks=False 

    for frame in [frames[0]]:
    #for frame in frames:    
        for sim in sims:
        #for sim in ['High_512']:
            beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]
            for param in params:
            #for param in [[1,1]]:
                at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                data_loc = at.get_data_location()
                for axis in axes:              
                #for axis in ['x']:

                    print "\n\n******* frame %s  %s  n0-%04d  %s *******\n\n"%(frame, beta, param[0], axis)
                    sys.stdout.flush()

                    field_horizontal = {'x':'By','y':'Bz','z':'Bx'}[axis]
                    field_vertical   = {'x':'Bz','y':'Bx','z':'By'}[axis]  
                    
                    # Get frbs
                    aQ = at.get_fits_array('Q%s_n0-%04d_p-%d'%(axis,param[0],param[1]), axis)
                    aU = at.get_fits_array('U%s_n0-%04d_p-%d'%(axis,param[0],param[1]), axis)
                    aN2 = at.get_fits_array('N2%s_n0-%04d_p-%d'%(axis,param[0],param[1]), axis)                            
                    aN = at.get_fits_array('density',axis)
                    aB_hor = at.get_fits_array(field_horizontal, axis)
                    aB_vert= at.get_fits_array(field_vertical, axis)                 
                    dg = at.get_dendrogram('density', axis)

                    # These lists below from pFrac_leaves to mass_leaves are sets of data in each disk disks for each leaf
                    # These get reset for each new axis, beta, param, and frame
                    # they have the form
                    # 
                    #            data = [ [leaf_0], [leaf_1], ..., [leaf_N] ]
                    #  where [leaf_N] = [ disk_sm, disk_1, ... , disk_N ]
                    #
                    #  where disk_N is a value in disk_N, for instance the average magnetic field angle in disk_N                    
                    pFrac_leaves = []              # polarization fraction in each disk
                    P_leaves_theta = []            # polarization angle in each disk
                    B_leaves_theta = []            # magnetic field angle in each disk
                    P_leaves_delta_theta = []      # |P_theta_large_disk - P_theta_small_disk| for each large disk
                    B_leaves_delta_theta = []      # |B_theta_large_disk - B_theta_small_disk| for each large disk
                    adjusted_P_leaves_theta = []   # polarization angle in disk minus mean field angle for the whole box 
                    adjusted_B_leaves_theta = []   # magnetic field angle in disk minus mean field angle for the whole box 
                    mass_leaves = []               # mass in each disk

                    r_lg_pc_list = []     # list of disk radii for each each disk
                    leafidx = []          # list of leaf indexes in the density dendrogram 
                    leaf_len = []             # list of length of leaf in pixels


                    BOXSIZE_PC = 4.6                 # Simulation is 4.6 parsecs on one side
                    PC_TO_PIXEL = res/BOXSIZE_PC     # [pixels/parsec]
                    p0 = 0.1     # Intrinsic normalization factor for polarization fraction
                                 # determined by the effective polarization cross section
                                 # for dust grains. p0 ~0.1 from observations. 
    
                    # Diameter of 0.1 pc (large) and 0.01 pc (small) disks for 8192 res
                    # or 1 and 0.1 for 512 res 
                    if res == 8192:
                        r_sm_pc = 0.005
                        r_lg_pc_max = 0.05            
                    elif res == 512:
                        r_sm_pc = 0.05
                        r_lg_pc_max = 0.5
                    else:
                        print '**** res should be 8192 or 512 for 512^3 datasets ****\nReturning Now.'
                        return -1
    
                    # Convert disk sized in pc to integer sizes in pixels
                    r_sm = int(round(r_sm_pc * PC_TO_PIXEL))
                    r_lg_max = int(round(r_lg_pc_max * PC_TO_PIXEL))
    
                    # Make list of large disk radii
                    if theta_vs_r_whole_box:    
                        r_lg_pc_list = np.array([r_sm_pc, 0.05, 0.1, 0.25, 0.5, 1., 2.])
                    else:
                        num_lg_disks = 25     # Number of disks from small scale to large scale
                        r_lg_pc_list = np.linspace(r_sm_pc, r_lg_pc_max, num_lg_disks)
                    
                    r_lg_list = (r_lg_pc_list*PC_TO_PIXEL).round().astype(int)           
                    BOX_CENTER = (res/2, res/2)    # center pixel of the projected simulation box
    
                    Btheta_boxAvg = np.mean(np.arctan2(aB_vert,aB_hor))*180/np.pi
    
                    k=0    # leaf counter
                    for leaf in dg.leaves: 
                        # These lists get reset for each leaf, and appended to for each large disk
                        pol_frac = []      # Polarization fraction large disk
                        Pangle_sm = []     #
                        Pangle_lg = []     #   
                        Bangle_sm = []     # 
                        Bangle_lg = []     # 
                        mass_in_disk = []  #

                        leafCenter = leaf.get_peak()[0]      # Pixel position of center of leaf in dendrogram [y,x]
    
                        # Roll the frbs so that the center of the leaf is in the center of the frb
                        # roll up or down (shift along y axis)
                        leaf_Q = np.roll(aQ, BOX_CENTER[0] - leafCenter[0], axis=0)
                        leaf_U = np.roll(aU, BOX_CENTER[0] - leafCenter[0], axis=0) 
                        leaf_N = np.roll(aN, BOX_CENTER[0] - leafCenter[0], axis=0)
                        leaf_N2 = np.roll(aN2, BOX_CENTER[0] - leafCenter[0], axis=0)  
                        leaf_B_hor = np.roll(aB_hor, BOX_CENTER[0] - leafCenter[0], axis=0)
                        leaf_B_vert = np.roll(aB_vert, BOX_CENTER[0] - leafCenter[0], axis=0)  
    
                        # roll left or right (shift along x axis)
                        leaf_Q = np.roll(leaf_Q, BOX_CENTER[1] - leafCenter[1], axis=1)
                        leaf_U = np.roll(leaf_U, BOX_CENTER[1] - leafCenter[1], axis=1) 
                        leaf_N = np.roll(leaf_N, BOX_CENTER[1] - leafCenter[1], axis=1)
                        leaf_N2 = np.roll(leaf_N2, BOX_CENTER[1] - leafCenter[1], axis=1)  
                        leaf_B_hor = np.roll(leaf_B_hor, BOX_CENTER[1] - leafCenter[1], axis=1)
                        leaf_B_vert = np.roll(leaf_B_vert, BOX_CENTER[1] - leafCenter[1], axis=1)              

                        leafidx.append(leaf.idx)             # Add leaf index to list of leaf indices  
    
                        # Find approximate size of leaf
                        leaf_len.append(leaf.indices()[0].max() - leaf.indices()[0].min())
                        leaf_len.append(leaf.indices()[1].max() - leaf.indices()[1].min())
        
                        n = 1     # large disk counter

                        for r_lg in r_lg_list:
                            # Cut out a square of data around leaf peak. Square size = 2*r_lg
                            cutoutSize = 2*r_lg
                            cutoutCenter = [cutoutSize/2,cutoutSize/2]            
                            ymin = BOX_CENTER[0] - cutoutSize/2     # Edges of square cutout
                            ymax = BOX_CENTER[0] + cutoutSize/2     #
                            xmin = BOX_CENTER[1] - cutoutSize/2     #
                            xmax = BOX_CENTER[1] + cutoutSize/2     #
    
                            y, x = np.mgrid[0:cutoutSize:1, 0:cutoutSize:1]    # Greate integer grid same resolution as the box cutout
    
                            # Calculate angle difference between small disk and each large disk for each leaf
                            # Make masks for the disks
                            mask_sm = (x - cutoutCenter[1])**2 + (y - cutoutCenter[0])**2 < r_sm**2
                            mask_lg = (x - cutoutCenter[1])**2 + (y - cutoutCenter[0])**2 < r_lg**2
                            
                            if total_masks and not theta_vs_r:
                                total_mask_lg = np.logical_or(mask_lg, total_mask_lg)    # Make one mask for all leaves 
                                total_mask_sm = np.logical_or(mask_sm, total_mask_sm)   
                                
                            # Cut out a square of data around leaf peak.    
                            Q = leaf_Q[ymin:ymax, xmin:xmax]
                            U = leaf_U[ymin:ymax, xmin:xmax] 
                            N = leaf_N[ymin:ymax, xmin:xmax]
                            N2 = leaf_N2[ymin:ymax, xmin:xmax]  
                            B_hor = leaf_B_hor[ymin:ymax, xmin:xmax]
                            B_vert = leaf_B_vert[ymin:ymax, xmin:xmax]                
                    
                            # Mask and average data over small and large scale disks
                            Q_avg_sm = np.mean(mask_sm * Q)
                            Q_avg_lg = np.mean(mask_lg * Q)
                            U_avg_sm = np.mean(mask_sm * U)
                            U_avg_lg = np.mean(mask_lg * U)
                            N_avg_sm = np.mean(mask_sm * N)
                            N_avg_lg = np.mean(mask_lg * N)
                            N2_avg_sm = np.mean(mask_sm * N2)
                            N2_avg_lg = np.mean(mask_lg * N2)  
                            B_hor_avg_sm = np.mean(mask_sm * B_hor)
                            B_hor_avg_lg = np.mean(mask_lg * B_hor)
                            B_vert_avg_sm = np.mean(mask_sm * B_vert)
                            B_vert_avg_lg = np.mean(mask_lg * B_vert)
                            mass_sm = np.sum(mask_sm * N)                  
                            mass_lg = np.sum(mask_lg * N)
    
                            # Calculate polarization angle and B field angle for small and large scale
                            Pangle_sm.append(0.5*np.arctan2(U_avg_sm,Q_avg_sm)*180/np.pi)
                            Pangle_lg.append(0.5*np.arctan2(U_avg_lg,Q_avg_lg)*180/np.pi)
                            Bangle_sm.append(np.arctan2(B_vert_avg_sm, B_hor_avg_sm)*180/np.pi)
                            Bangle_lg.append(np.arctan2(B_vert_avg_lg, B_hor_avg_lg)*180/np.pi)
                            mass_in_disk.append(mass_lg) 
                            pol_frac.append( p0 * np.sqrt(Q_avg_lg**2 + U_avg_lg**2)/(N_avg_lg - p0*N2_avg_lg) )                                          
                                
                            n+=1    # Increment large disk counter
                            ## End for r_lg in r_lg_list ##
    
                        # Calculate angle differences between small disk and all large disks  
                        # and add data to lists for a single leaf  
                        Pangle_sm = np.array(Pangle_sm)
                        Pangle_lg = np.array(Pangle_lg)
                        Bangle_sm = np.array(Bangle_sm)
                        Bangle_lg = np.array(Bangle_lg)                            
                        Pangle_diff = np.abs(Pangle_lg - Pangle_sm)
                        Bangle_diff = np.abs(Bangle_lg - Bangle_sm)                            
                        P_leaves_delta_theta.append(Pangle_diff)
                        B_leaves_delta_theta.append(Bangle_diff)
                        P_leaves_theta.append(Pangle_lg)
                        B_leaves_theta.append(Bangle_lg)
                        pFrac_leaves.append(pol_frac)
                        mass_leaves.append(mass_in_disk)
    
                        k+=1    # Increment leaf counter
                        ## End for leaf in leaves: ##
    
                    avg_leaf_size = np.mean(np.array(leaf_len))
                    #print "\nApproximate leaf size:", np.round(avg_leaf_size,3), "pixels"
                    #print "Approximate leaf size:", np.round(avg_leaf_size*1/PC_TO_PIXEL,3), 'pc'
                    
                    if total_masks and not theta_vs_r:
                        density_lg = total_mask_lg * N
                        density_sm = total_mask_sm * N
                    
                    # convert lists to arrays
                    P_leaves_delta_theta = np.array(P_leaves_delta_theta)
                    B_leaves_delta_theta = np.array(B_leaves_delta_theta)
                    P_leaves_theta = np.array(P_leaves_theta)
                    B_leaves_theta = np.array(B_leaves_theta)
                    pFrac_leaves = np.array(pFrac_leaves)
                    mass_leaves = np.array(mass_leaves)
                    adjusted_P_leaves_theta = np.abs(P_leaves_theta - Btheta_boxAvg)
                    adjusted_B_leaves_theta = np.abs(B_leaves_theta - Btheta_boxAvg)
         
                    # Remove degeneracy from delta_theta        
                    # If delta_theta is > 90 deg, use 180 minus delta_theta                            
                    P_leaves_delta_theta[P_leaves_delta_theta > 90] = 180 - P_leaves_delta_theta[P_leaves_delta_theta > 90]  
                    B_leaves_delta_theta[B_leaves_delta_theta > 180] = 360 - B_leaves_delta_theta[B_leaves_delta_theta > 180]
                    B_leaves_delta_theta[B_leaves_delta_theta > 90] = 180 - B_leaves_delta_theta[B_leaves_delta_theta > 90]                        
                    adjusted_P_leaves_theta[adjusted_P_leaves_theta > 90] = 180 - adjusted_P_leaves_theta[adjusted_P_leaves_theta > 90]  
                    adjusted_B_leaves_theta[adjusted_B_leaves_theta > 180] = 360 - adjusted_B_leaves_theta[adjusted_B_leaves_theta > 180]
                    adjusted_B_leaves_theta[adjusted_B_leaves_theta > 90] = 180 - adjusted_B_leaves_theta[adjusted_B_leaves_theta > 90]

                    # Get data for each leaf's small disk 
                    pol_frac_sm = np.array([item[0] for item in pFrac_leaves])
                    P_theta_sm = np.array([item[0] for item in P_leaves_theta])
                    B_theta_sm = np.array([item[0] for item in B_leaves_theta])
                    mass_in_sm_disk = np.array([item[0] for item in mass_leaves])

                    # Find maximum polarization fraction out of all axes, beta, n0, and frame
                    if pol_frac_sm.max() > pFrac_max:
                        pFrac_max = pol_frac_sm.max()
                    
                    # Pickle data 
                    if theta_vs_r_whole_box:
                        bdump(pFrac_max, "%s/P_frac_max_whole_box_frame%04d.data" %(data_loc, frame))
                        # Each of the below arrays has the form 
                        #            data = [ [leaf_0], [leaf_1], ..., [leaf_N] ]
                        #  where [leaf_N] = [ disk_sm, disk_1, ... , disk_N ]
                        bdump(pFrac_leaves, "%s/P_frac_leaves_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))                                         
                        bdump(P_leaves_theta, "%s/P_leaves_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                        bdump(B_leaves_theta, "%s/B_leaves_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))                         
                        bdump(P_leaves_delta_theta, "%s/P_leaves_delta_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                        bdump(B_leaves_delta_theta, "%s/B_leaves_delta_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                        bdump(adjusted_P_leaves_theta, "%s/adjusted_P_leaves_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                        bdump(adjusted_B_leaves_theta, "%s/adjusted_B_leaves_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta)) 
                        bdump(r_lg_pc_list, "%s/r_lg_pc_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                        bdump(mass_leaves, "%s/mass_leaves_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))                        
                    else:
                       bdump(pFrac_max, "%s/P_frac_max_frame%04d.data" %(data_loc, frame))                                                
                       bdump(pFrac_leaves, "%s/P_frac_leaves_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))                                         
                       bdump(P_leaves_theta, "%s/P_leaves_theta_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                       bdump(B_leaves_theta, "%s/B_leaves_theta_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))                         
                       bdump(P_leaves_delta_theta, "%s/P_leaves_delta_theta_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                       bdump(B_leaves_delta_theta, "%s/B_leaves_delta_theta_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                       bdump(adjusted_P_leaves_theta, "%s/adjusted_P_leaves_theta_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                       bdump(adjusted_B_leaves_theta, "%s/adjusted_B_leaves_theta_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta)) 
                       bdump(r_lg_pc_list, "%s/r_lg_pc_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                       bdump(mass_leaves, "%s/mass_leaves_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))

                    # Pickle just small disk data for each leaf
                    # data = [ leaf_0_disk_sm, leaf_1_disk_sm, ... , leaf_N_disk_sm]
                    bdump(leafidx, "%s/leafidx_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))                                            
                    bdump(pol_frac_sm, "%s/P_frac_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))                                                         
                    bdump(P_theta_sm, "%s/P_theta_sm_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                    bdump(B_theta_sm, "%s/Btheta_sm_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                    bdump(mass_in_sm_disk, "%s/mass_in_sm_disk_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                                 
#---- end analyze_dense_regions ----#

# Make P_delta_theta vs P_frac plots
def make_P_frac_plots():
    """
    make scatter plot of polarization delta_theta vs polarization fraction
    """
    #for frame in [frames[2]]:
    for frame in frames:
        #param has the form [n0,p]
        #for param in [[1,1]]:    
        for param in params:
            #for axis in ['x']:        
            for axis in axes: 
                # Lists of data element 0, 1, 2 are B02, B2, B20 respectively                     
                #for sim in ['High_512']:
                for sim in sims:
                    beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]               
                
                    at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                    plot_loc_Pfrac = at.get_plot_location_P_frac()
                    data_loc = at.get_data_location()
                    
                    # Load pickled data
                    pFrac = bload("%s/P_frac_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                    P_leaves_delta_theta = bload("%s/P_leaves_delta_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                    pFrac_max = bload("%s/P_frac_max_whole_box_frame%04d.data" %(data_loc, frame))

                    for i in range(1,len(r_disk)):
                        P_delta_theta = np.array([item[i] for item in P_leaves_delta_theta])

                        #fit = np.polyfit(pFrac, P_delta_theta, deg=1)                
             
                        # Make plot of Polarization Angle vs Polarization Fraction
                        print "plotting angle vs Pfrac", axis, beta, "n0 =", param[0], 'p =', param[1], r_disk[i], 'pc disk' 
                        plt.clf()
                        plt.scatter(pFrac, P_delta_theta)
                        #plt.plot(pFrac,fit[0]*pFrac + fit[1], color='red')
                        plt.title(r'$\Delta\theta$'+' vs Polarization Fraction %s %s n0 = %04d, p = %d\n%s pc disk'%(axis, beta, param[0], param[1], r_disk[i]))
                        plt.xlabel(r'$\bar{P}_{frac}$')
                        plt.ylabel(r'$|\theta_{lg} - \theta_{sm}|$', size=15)
                        #plt.axis([0,pFrac_max,0,90])
                        plt.axis([0,0.072,0,90])
                        plt.grid()
                        plt.savefig("%s/PAngle_vs_Pfraction_frame%04d_%s_n0-%04d_p-%d_%s_disk-%s" %(plot_loc_Pfrac, frame, beta, param[0], param[1], axis, r_disk_str[i]))  

# Make plots while varying field strength
def make_field_strength_plots():
    """
    plot histograms of polarization delta_theta and magnetic field delta_theta.
    put three histograms on one plot while varying the initial mean magnetic field
    strength. Plasma beta = 0.2, 2, and 20 (B02, B2, B20 respectively)
    """  
    #for frame in [frames[2]]:  
    for frame in frames:
        #param has the form [n0,p]
        #for param in [[1,1]]:    
        for param in params:
            #for axis in ['x']:        
            for axis in axes: 
                for i in range(len(r_disk)):
                    # Lists of data element 0, 1, 2 are B02, B2, B20 respectively
                    P_delta_theta = []
                    B_delta_theta = []
                    adjusted_Btheta = []
                    adjusted_P_theta = []  

                    weights = []
                     
                    for sim in sims:
                        beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]               
                    
                        at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                        plot_loc2 = at.get_plot_location_field_strength()
                        data_loc = at.get_data_location()
      
                        # Load pickled data
                        P_leaves_delta_theta = bload("%s/P_leaves_delta_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                        B_leaves_delta_theta = bload("%s/B_leaves_delta_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                        adjusted_P_leaves_theta = bload("%s/adjusted_P_leaves_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                        adjusted_B_leaves_theta = bload("%s/adjusted_B_leaves_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta)) 
                            
                        # Extract data in disk disk_i for each leaf    
                        P_delta_theta.append([item[i] for item in P_leaves_delta_theta])
                        B_delta_theta.append([item[i] for item in B_leaves_delta_theta])
                        adjusted_P_theta.append([item[i] for item in adjusted_P_leaves_theta])                        
                        adjusted_Btheta.append([item[i] for item in adjusted_B_leaves_theta])
                        
                        ## Weight each histogram by number of leaves
                        ##data = [item[i] for item in P_leaves_delta_theta]
                        #weights.append(list(np.ones_like(data)*len(data)))
                    '''
                    print "\nweights", weights
                    weights_max = max(max(weights))
                    print "\nweights max", weights_max
                    weights[:] = [weights_max/x for x in weights]
                    print "\nweights", weights
                    '''
                    # Histogram of polarization angle difference
                    color = ['red', 'green', 'blue']
                    plt.clf()
                    plt.title(r'$\Delta\theta$'+' %s n0 = %04d, p = %d frame %d %s pc disk'%(axis, param[0], param[1], frame, r_disk[i]))
                    plt.xlabel(r'$|\theta_{lg} - \theta_{sm}|$', size=15)
                    plt.xlim(0,90)
                    #plt.hist(P_delta_theta, histtype='step', color=color, label=['B02','B2','B20'], bins=bins, weights=weights, normed=True)
                    plt.hist(P_delta_theta, histtype='step', color=color, label=['B02','B2','B20'], bins=bins)
                    plt.legend()
                    plt.savefig("%s/P_angle_diff_histogram_frame%04d_n0-%04d_p-%d_%s_disk-%s_field_strength" %(plot_loc2, frame, param[0], param[1], axis, r_disk_str[i]))
                    
                    # Histogram of B field angle difference
                    plt.clf()
                    plt.title(r'$\Delta\theta$'+' Magnetic Field %s axis frame %d %s pc disk'%(axis, frame, r_disk[i]))
                    plt.xlabel(r'$|\theta_{lg} - \theta_{sm}|$', size=15)
                    plt.xlim(0,90)
                    #plt.hist(B_delta_theta, histtype='step', color=color, label=['B02','B2','B20'], bins=bins, weights=weights, normed=True)
                    plt.hist(B_delta_theta, histtype='step', color=color, label=['B02','B2','B20'], bins=bins)
                    plt.legend()
                    plt.savefig("%s/Bfield_angle_diff_histogram_frame%04d_%s_disk-%s_field_strength" %(plot_loc2, frame, axis, r_disk_str[i])) 
        
                    # Histogram of B field angle large minus mean B field angle
                    plt.clf()
                    plt.title('B Field Angle of Large Disk Minus Mean B Field Angle of Box\n %s axis frame %d %s pc disk'%(axis, frame, r_disk[i]))
                    plt.xlabel(r'$|\theta_{lg} - \theta_{BoxAvg}|$', size=15)
                    plt.xlim(0,90)
                    #plt.hist(adjusted_Btheta, histtype='step', color=color, label=['B02','B2','B20'], bins=bins, weights=weights, normed=True)
                    plt.hist(adjusted_Btheta, histtype='step', color=color, label=['B02','B2','B20'], bins=bins)
                    plt.legend()
                    plt.savefig("%s/adjusted_Bfield_angle_histogram_frame%04d_%s_disk-%s_field_strength" %(plot_loc2, frame, axis, r_disk_str[i]))
        
                    # Histogram of Polarization angle large minus mean B field angle
                    plt.clf()
                    plt.title('Polarization Angle of Large Disk Minus Mean B Field Angle of Box\n%s n0 = %04d, p = %d frame %d %s pc disk'%(axis, param[0], param[1], frame, r_disk[i]))
                    plt.xlabel(r'$|\theta_{lg} - \theta_{BoxAvg}|$', size=15)
                    plt.xlim(0,90)
                    #plt.hist(adjusted_P_theta, histtype='step', color=color, label=['B02','B2','B20'], bins=bins, weights=weights, normed=True)
                    plt.hist(adjusted_P_theta, histtype='step', color=color, label=['B02','B2','B20'], bins=bins)
                    plt.legend()
                    plt.savefig("%s/adjusted_P_angle_histogram_frame%04d_n0-%04d_p-%d_%s_disk-%s_field_strength" %(plot_loc2, frame, param[0], param[1], axis, r_disk_str[i]))
                    

# Make plots with variying viewing angle       
def make_viewing_angle_plots():
    """
    plot histograms of polarization delta_theta and magnetic field delta_theta.
    put three histograms on one plot while varying the the viewing angle x, y, and z.
    """ 
    #for frame in [frames[2]]:
    for frame in frames:
        #for param in [[1,1]]:
        for param in params:  
            for sim in sims:
                for i in range(len(r_disk)):

                    # Lists of data element 0, 1, 2 are for x, y, z respectively
                    P_delta_theta = []
                    B_delta_theta = []
                    adjusted_Btheta_sm = []
                    adjusted_Btheta = []
                    adjusted_P_theta_sm = []
                    adjusted_P_theta = []
                    weights = []
    
                    beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]  
                    at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                    plot_loc2 = at.get_plot_location_viewing_angle()
                    data_loc = at.get_data_location()      
    
                    for axis in axes:                                
                        # Load pickled data
                        P_leaves_delta_theta = bload("%s/P_leaves_delta_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                        B_leaves_delta_theta = bload("%s/B_leaves_delta_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                        adjusted_P_leaves_theta = bload("%s/adjusted_P_leaves_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                        adjusted_B_leaves_theta = bload("%s/adjusted_B_leaves_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta)) 
                            
                        P_delta_theta.append([item[i] for item in P_leaves_delta_theta])
                        B_delta_theta.append([item[i] for item in B_leaves_delta_theta])
                        adjusted_P_theta.append([item[i] for item in adjusted_P_leaves_theta])                        
                        adjusted_Btheta.append([item[i] for item in adjusted_B_leaves_theta])

                        # Weight each histogram by number of leaves
                        #data = [item[i] for item in P_leaves_delta_theta]
                        #weights.append(list(np.ones_like(data)*len(data)))                        
    
                    # Histogram of polarization angle difference
                    color = ['darkred', 'steelblue', 'olivedrab']            
                    plt.clf()
                    plt.title(r'$\Delta\theta$'+' %s n0 = %04d, p = %d frame %d\n%s pc disk'%(beta, param[0], param[1], frame, r_disk[i]))
                    plt.xlabel(r'$|\theta_{lg} - \theta_{sm}|$', size=15)
                    plt.xlim(0,90)
                    #plt.hist(P_delta_theta, histtype='step', color=color, label=axes, bins=bins, weights=weights, normed=True)
                    plt.hist(P_delta_theta, histtype='step', color=color, label=axes, bins=bins)
                    plt.legend()
                    plt.savefig("%s/P_angle_diff_histogram_frame%04d_n0-%04d_p-%d_%s_disk-%s_viewing_angle" %(plot_loc2, frame, param[0], param[1], beta, r_disk_str[i]))

                    # Histogram of B field angle difference
                    plt.clf()
                    plt.title(r'$\Delta\theta$'+' Magnetic Field %s frame %d\n%s pc disk'%(beta, frame, r_disk[i]))
                    plt.xlabel(r'$|\theta_{lg} - \theta_{sm}|$', size=15)
                    plt.xlim(0,90)
                    #plt.hist(B_delta_theta, histtype='step', color=color, label=axes, bins=bins, weights=weights, normed=True)
                    plt.hist(B_delta_theta, histtype='step', color=color, label=axes, bins=bins)
                    plt.legend()
                    plt.savefig("%s/Bfield_angle_diff_histogram_frame%04d_%s_disk-%s_viewing_angle" %(plot_loc2, frame, beta, r_disk_str[i])) 
    
                    # Histogram of B field angle large minus mean B field angle
                    plt.clf()
                    plt.title('B Field Angle of Large Disk Minus Mean B Field Angle of Box\n%s frame %d %s pc disk'%(beta, frame, r_disk[i]))
                    plt.xlabel(r'$|\theta_{lg} - \theta_{BoxAvg}|$', size=15)
                    plt.xlim(0,90)
                    #plt.hist(adjusted_Btheta, histtype='step', color=color, label=axes, bins=bins, weights=weights, normed=True)
                    plt.hist(adjusted_Btheta, histtype='step', color=color, label=axes, bins=bins)
                    plt.legend()
                    plt.savefig("%s/adjusted_Bfield_angle_histogram_frame%04d_%s_disk-%s_viewing_angle" %(plot_loc2, frame, beta, r_disk_str[i]))
    
                    # Histogram of polarization angle large minus mean B field angle
                    plt.clf()
                    plt.title('Polarization Angle of Large Disk Minus Mean B Field Angle of Box\n%s n0 = %04d, p = %d frame %d %s pc disk'%(beta, param[0], param[1], frame, r_disk[i]))
                    plt.xlabel(r'$|\theta_{lg} - \theta_{BoxAvg}|$', size=15)
                    plt.xlim(0,90)
                    #plt.hist(adjusted_P_theta, histtype='step', color=color, label=axes, bins=bins, weights=weights, normed=True)
                    plt.hist(adjusted_P_theta, histtype='step', color=color, label=axes, bins=bins)
                    plt.legend()
                    plt.savefig("%s/adjusted_P_angle_histogram_frame%04d_n0-%04d_p-%d_%s_disk-%s_viewing_angle" %(plot_loc2, frame, param[0], param[1], beta, r_disk_str[i]))

def make_time_dependence_plots():
    """
    plot histograms of polarization delta_theta and magnetic field delta_theta.
    put three histograms on one plot while varying the frame. frames 10 30 and 50.
    """
    #for frame in [50]:
    for frame in frames:
        #for sim in [sims[0]]:
        for sim in sims:
            beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]
            for axis in axes:
                for i in range(len(r_disk)):
                    for param in [[1,1]]:
                    #for param in params:
                
                        P_delta_theta = []
                        B_delta_theta = []
                        adjusted_Btheta_sm = []
                        adjusted_Btheta = []
                        adjusted_P_theta_sm = []
                        adjusted_P_theta = []
                        weights = []
    
        
                        for frame in frames:
                                               
                            at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                            data_loc = at.get_data_location()
                            plot_loc = at.get_plot_location_time_dependence()
        
                            P_leaves_delta_theta = bload("%s/P_leaves_delta_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                            B_leaves_delta_theta = bload("%s/B_leaves_delta_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                            adjusted_P_leaves_theta = bload("%s/adjusted_P_leaves_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                            adjusted_B_leaves_theta = bload("%s/adjusted_B_leaves_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta)) 
                                
                            P_delta_theta.append([item[i] for item in P_leaves_delta_theta])
                            B_delta_theta.append([item[i] for item in B_leaves_delta_theta])
                            adjusted_P_theta.append([item[i] for item in adjusted_P_leaves_theta])                        
                            adjusted_Btheta.append([item[i] for item in adjusted_B_leaves_theta])

                            #data = [item[i] for item in P_leaves_delta_theta]
                            #weights.append(list(np.ones_like(data)*len(data)))
                         
                        # Histogram of polarization angle difference
                        color = ['gold', 'firebrick', 'indigo']            
                        plt.clf()
                        plt.title(r'$\Delta\theta$'+' %s %s n0 = %04d, p = %d %s pc disk'%(axis, beta, param[0], param[1], r_disk[i]))
                        plt.xlabel(r'$|\theta_{lg} - \theta_{sm}|$', size=15)
                        plt.xlim(0,90)
                        #plt.hist(P_delta_theta, histtype='step', color=color, label=["frame %04d"%frames[0], "frame %04d"%frames[1], "frame %04d"%frames[2]], bins=bins, weights=weights, normed=True)
                        plt.hist(P_delta_theta, histtype='step', color=color, label=["frame %04d"%frames[0], "frame %04d"%frames[1], "frame %04d"%frames[2]], bins=bins)
                        plt.legend()
                        plt.savefig("%s/P_angle_diff_histogram_%s_%s_n0-%04d_p-%d_disk-%s_time_dependence" %(plot_loc, axis, beta, param[0], param[1], r_disk_str[i]))
            
                        # Histogram of polarization angle large minus mean B field angle
                        plt.clf()
                        plt.title('Polarization Angle of Large Disk Minus Mean B Field Angle of Box\n%s %s n0 = %04d, p = %d %s pc disk'%(axis, beta, param[0], param[1], r_disk[i]))
                        plt.xlabel(r'$|\theta_{lg} - \theta_{BoxAvg}|$', size=15)
                        plt.xlim(0,90)
                        #plt.hist(adjusted_P_theta, histtype='step', color=color, label=["frame %04d"%frames[0], "frame %04d"%frames[1], "frame %04d"%frames[2]], bins=bins, weights=weights, normed=True)
                        plt.hist(adjusted_P_theta, histtype='step', color=color, label=["frame %04d"%frames[0], "frame %04d"%frames[1], "frame %04d"%frames[2]], bins=bins)
                        plt.legend()
                        plt.savefig("%s/adjusted_P_angle_histogram_%s_%s_n0-%04d_p-%d_disk-%s_time_dependence" %(plot_loc, axis, beta, param[0], param[1], r_disk_str[i]))
                        
                    # Histogram of B field angle difference
                    plt.clf()
                    plt.title(r'$\Delta\theta$'+' Magnetic Field %s %s %s pc disk'%(axis, beta, r_disk[i]))
                    plt.xlabel(r'$|\theta_{lg} - \theta_{sm}|$', size=15)
                    plt.xlim(0,90)
                    #plt.hist(B_delta_theta, histtype='step', color=color, label=["frame %04d"%frames[0], "frame %04d"%frames[1], "frame %04d"%frames[2]], bins=bins, weights=weights, normed=True)
                    plt.hist(B_delta_theta, histtype='step', color=color, label=["frame %04d"%frames[0], "frame %04d"%frames[1], "frame %04d"%frames[2]], bins=bins)
                    plt.legend()
                    plt.savefig("%s/Bfield_angle_diff_histogram_%s_%s_disk-%s_time_dependence" %(plot_loc, axis, beta, r_disk_str[i]))
        
                    # Histogram of B field angle large minus mean B field angle
                    plt.clf()
                    plt.title('B Field Angle of Large Disk Minus Mean\nB Field Angle of Box %s %s %s pc disk'%(axis, beta, r_disk[i]))
                    plt.xlabel(r'$|\theta_{lg} - \theta_{BoxAvg}|$', size=15)
                    plt.xlim(0,90)
                    #plt.hist(adjusted_Btheta, histtype='step', color=color, label=["frame %04d"%frames[0], "frame %04d"%frames[1], "frame %04d"%frames[2]], bins=bins, weights=weights, normed=True)
                    plt.hist(adjusted_Btheta, histtype='step', color=color, label=["frame %04d"%frames[0], "frame %04d"%frames[1], "frame %04d"%frames[2]], bins=bins)
                    plt.legend()
                    plt.savefig("%s/adjusted_Bfield_angle_histogram_%s_%s_disk-%s_time_dependence" %(plot_loc, axis, beta, r_disk_str[i]))

# Make plots with variying threshold density      
def make_threshold_density_plots():
    """
    plot histograms of polarization delta_theta and magnetic field delta_theta.
    put three histograms on one plot while varying the the cutoff density above which
    there is depolarization. Cutoff densities n0 = 19, 39 and 1945.
    """
    params = [[19,0], [39,0], [1945,0]]
    #for frame in [50]:
    for frame in frames:
        for axis in axes:   
            for sim in sims:
                # List of data, element 0, 1, 2 are for n0 = 39, 19, 1945 respectively
                beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim] 
                
                for i in range(len(r_disk)):     
                    P_delta_theta = []
                    for param in params:  
        
                        at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                        plot_loc = at.get_plot_location_threshold_density()
                        data_loc = at.get_data_location()    
                                                               
                        # Load pickled data
                        P_leaves_delta_theta = bload("%s/P_leaves_delta_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                        P_delta_theta.append([item[i] for item in P_leaves_delta_theta])
    
                    # Histogram of polarization angle difference
                    plt.clf()
                    plt.title(r'$\Delta\theta$'+' %s %s frame %d %s pc disk'%(axis, beta, frame, r_disk[i]))
                    plt.xlabel(r'$|\theta_{lg} - \theta_{sm}|$', size=15)
                    plt.xlim(0,90)
                    color = ['lightpink', 'mediumvioletred','darkslateblue']
                    plt.hist(P_delta_theta, histtype='step', color=color, label=['n0 = %04d'%params[0][0],
                                                                                 'n0 = %04d'%params[1][0],
                                                                                 'n0 = %04d'%params[2][0]], bins=bins)
                    plt.legend()
                    plt.savefig("%s/P_angle_diff_histogram_frame%04d_%s_%s_disk-%s_threshold_density" %(plot_loc, frame, beta, axis, r_disk_str[i]))
   

def make_P_delta_theta_vs_B_delta_theta_plots():
    """
    Compare polarization delta_theta and magnetic field delta theta by 
    plotting histograms of each on one plot and making a scatter plot of 
    P_delta_theta vs B_delta_theta.
    """
    #for frame in [50]:
    for frame in frames:
        #for sim in [sims[0]]:
        for sim in sims:
    
            beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]
    
            for param in [[1,1]]:
    
                at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                plot_loc = at.get_plot_location_P_theta_vs_B_theta()
                data_loc = at.get_data_location() 
    
                for axis in axes:                    
                    P_leaves_theta = bload("%s/P_leaves_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                    B_leaves_theta = bload("%s/B_leaves_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                    P_leaves_delta_theta = bload("%s/P_leaves_delta_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                    B_leaves_delta_theta = bload("%s/B_leaves_delta_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))

                    # convert angles in range [90,180] U [-90, -180] deg to angles in range [-90, 90] deg
                    B_leaves_theta[B_leaves_theta > 90] = B_leaves_theta[B_leaves_theta > 90] - 180
                    B_leaves_theta[B_leaves_theta < -90] = 180 + B_leaves_theta[B_leaves_theta < -90]

                    for i in range(len(r_disk)):
                        P_theta = [item[i] for item in P_leaves_theta]
                        B_theta = [item[i] for item in B_leaves_theta]
                        P_delta_theta = [item[i] for item in P_leaves_delta_theta]
                        B_delta_theta = [item[i] for item in B_leaves_delta_theta]
                        
                        data = [P_delta_theta, B_delta_theta]
         
                        # Histogram of angle differences
                        plt.clf()
                        plt.title(r'$\Delta\theta$'+'%s %s n0 = %04d, p = %d frame %d %s pc disk'%(axis, beta, param[0], param[1], frame, r_disk[i]))
                        plt.xlabel(r'$|\theta_{lg} - \theta_{sm}|$', size=15)
                        plt.xlim(0,90)
                        plt.hist(data, histtype='step', color=['blue','black'], label=['Polarization angle','Magnetic field angle'], bins=bins)
                        plt.legend()
                        plt.savefig("%s/P_and_B_angle_diff_histogram_frame%04d_n0-%04d_p-%d_%s_%s_disk-%s" %(plot_loc, frame, param[0], param[1], axis, beta, r_disk_str[i]))
        
                        # Scatter plot of angle differences
                        plt.clf()
                        plt.title(r'$\Delta\theta_{polarization}$ vs $\Delta\theta_{Magnetic}$'+'\n%s %s n0 = %04d, p = %d frame %d %s pc disk'%(axis, beta, param[0], param[1], frame, r_disk[i]))
                        plt.xlabel(r'$|\theta_{lg} - \theta_{sm}|'+' Magnetic Field$', size=15)
                        plt.ylabel(r'$|\theta_{lg} - \theta_{sm}|'+' Polarization$', size=15)              
                        plt.xlim(0,90)
                        plt.ylim(0,90)
                        plt.scatter(B_delta_theta, P_delta_theta)
                        plt.savefig("%s/P_and_B_angle_diff_scatter_frame%04d_n0-%04d_p-%d_%s_%s_disk-%s" %(plot_loc, frame, param[0], param[1], axis, beta, r_disk_str[i]))  
    
                        # Scatter plot of P and B angles
                        plt.clf()
                        plt.title(r'$\theta_{polarization}$ vs $\theta_{Magnetic}$'+'\n%s %s n0 = %04d, p = %d frame %d %s pc disk'%(axis, beta, param[0], param[1], frame, r_disk[i]))
                        plt.xlabel(r'$\theta$'+'  Magnetic Field', size=15)
                        plt.ylabel(r'$\theta$'+'  Polarization', size=15)              
                        plt.xlim(-95,95)
                        plt.ylim(-95,95)
                        plt.xticks(range(-90,91,30), range(-90,91,30))
                        plt.yticks(range(-90,91,30), range(-90,91,30))
                        plt.scatter(B_theta, P_theta)
                        plt.savefig("%s/P_and_B_theta_scatter_frame%04d_n0-%04d_p-%d_%s_%s_disk-%s" %(plot_loc, frame, param[0], param[1], axis, beta, r_disk_str[i]))

def make_theta_vs_r_plots():
    """
    plot polarization delta_theta and magnetic field delta theta as
    a function of radius of the larger disk.
    """
    for frame in [50]:
    #for frame in frames:
        # param has the form [n0,p]
        for param in [[1,1]]: 
    
            #for axis in ['x']:        
            for axis in axes:  
    
                #for sim in ['High_512']:
                for sim in sims:
    
                    beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]                
                    at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                    plot_loc_leaves = at.get_plot_location_leaves()
                    data_loc = at.get_data_location()
                    
                    if theta_vs_r_whole_box:
                        P_leaves_delta_theta = bload("%s/P_leaves_delta_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                        B_leaves_delta_theta = bload("%s/B_leaves_delta_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                        title = r'$\Delta\theta$'+' vs Radius of Large Disk Whole Box\n%s %s n0 = %04d, p = %d Leaf idx: %d frame %d'
                        save = "%s/leaf%03d_P_delta_theta_vs_r_disk_whole_box_%s_n0-%04d_p-%d_%s_frame%04d" 
                    else:                            
                        P_leaves_delta_theta = bload("%s/P_leaves_delta_theta_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                        B_leaves_delta_theta = bload("%s/B_leaves_delta_theta_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                        title = r'$\Delta\theta$'+' vs Radius of Large Disk\n%s %s n0 = %04d, p = %dLeaf idx: %d frame %d'
                        save = "%s/leaf%03d_P_delta_theta_vs_r_disk_%s_n0-%04d_p-%d_%s_frame%04d"

                    leafidx = bload("%s/leafidx_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                    
                    r_max = r_disk.max()
      
                    # Plot angle difference vs radius of large disk
                    # Polarization angle
                    print "plotting P_leaves_delta_theta %s %s %s %s"%(frame,param, axis, beta)
                    for i in range(len(P_leaves_delta_theta)):                  
                        plt.clf()
                        plt.title(title%(axis, beta, param[0], param[1], leafidx[i], frame))
                        plt.xlabel('Radius of Large Disk [pc]')
                        plt.ylabel(r'$|\theta_{lg} - \theta_{sm}|$', size=15)
                        plt.ylim(0,90)
                        plt.xlim(0,r_max)
                        plt.plot(r_disk, P_leaves_delta_theta[i])
                        plt.savefig(save%(plot_loc_leaves, leafidx[i], axis, param[0], param[1], beta, frame))
                    
                    # B field angle
                    '''
                    print "plotting B_leaves_delta_theta"
                    for i in range(len(B_leaves_delta_theta)):                 
                        plt.clf()
                        plt.title(r'$\Delta\theta$'+' Magnetic Field vs Radius of Large Disk %s %s \nLeaf idx: %d frame %d'%(axis, beta,leafidx[i], frame))
                        plt.xlabel('Radius of Large Disk [pc]')
                        plt.ylabel(r'$|\theta_{lg} - \theta_{sm}|$', size=15)
                        plt.ylim(0,90)
                        plt.xlim(0,r_max)                    
                        plt.plot(r_disk, B_leaves_delta_theta[i])
                        plt.savefig("%s/leaf%03d_Bfield_delta_theta_vs_r_disk_%s_%s_frame%04d" %(plot_loc_leaves, leafidx[i], axis, beta, frame))
                    '''


# This is for making density images with magnetic field overlays
def make_density_with_Bfield_plots(weight_pfrac=True):
    """
    Plot a image of the density of the whole simulation and an 
    image of each leaf. Also calcululate the polarization angle
    and polarization fraction and overlay the polars on top of
    the images. 
    """
    weight_pfrac=False
    for frame in [50]:
    #for frame in frames:
        for sim in sims:
        #for sim in ['High_512']:
            beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]
        
            # param has the form [n0,p]
            for param in [[1,1]]:        
    
                at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                plot_loc = at.get_plot_location_leaves()
        
                for axis in ['y']:
                #for axis in ['x', 'y', 'z']:
                    field_horizontal = {'x':'By','y':'Bz','z':'Bx'}[axis]
                    field_vertical   = {'x':'Bz','y':'Bx','z':'By'}[axis]                   
                    # Get .fits data
                    print "getting stokes_%s fits" %axis
                    Q = at.get_fits_array('Q%s_n0-%04d_p-%d'%(axis,param[0],param[1]), axis)
                    U = at.get_fits_array('U%s_n0-%04d_p-%d'%(axis,param[0],param[1]), axis)
                    N2 = at.get_fits_array('N2%s_n0-%04d_p-%d'%(axis,param[0],param[1]), axis)                            
                    N = at.get_fits_array('density',axis)    
                    B_hor = at.get_fits_array(field_horizontal, axis)
                    B_vert= at.get_fits_array(field_vertical, axis)

                    #dg = at.get_dendrogram('density', axis)
        
                    # Diameter of 0.1 pc (large) and 0.01 pc (small) disks for 8192 res
                    # or 1 and 0.1 for 512 res 
                    if res == 8192:
                        r_sm_pc = 0.005
                        r_lg_pc_max = 0.05            
                    elif res == 512:
                        r_sm_pc = 0.05
                        r_lg_pc_max = 0.5
                    else:
                        print '**** res should be 8192 or 512 for 512^3 datasets ****\nReturning Now.'
                        return -1
                    
                    # Convert disk sized in pc to integer sizes in pixels
                    BOXSIZE_PC = 4.6                 # Simulation is 4.6 parsecs on one side
                    PC_TO_PIXEL = res/BOXSIZE_PC     # [pixels/parsec]                    
                    r_sm = int(round(r_sm_pc * PC_TO_PIXEL))
                    r_lg_max = int(round(r_lg_pc_max * PC_TO_PIXEL))   
        
                    # Polarization angle
                    print "\ncalculating theta_stokes_%s" %axis
                    theta_stokes = 0.5*np.arctan2(U,Q)*180/np.pi
                    B_theta = np.arctan2(B_vert, B_hor)*180/np.pi 
                    #B_theta[B_theta > 90] = B_theta[B_theta > 90] - 180
                    #B_theta[B_theta < -90] = 180 + B_theta[B_theta < -90]                    
                    
                    # Calculate polarization fraction
                    print "calculating polarization_frac_%s" %axis
                    p0 = 0.1
                    polarization_frac = p0*np.sqrt(Q**2 + U**2)/(N - p0*N2)
                    
                    # Make vector field of polarization angles with length given by polarization fraction
                    SPACING = res/40   # space between each vector (40 vectors per side)
                
                    # grid of positions where to plot vectors
                    Y, X = np.mgrid[0:res:SPACING, 0:res:SPACING]
                    
                    if weight_pfrac:
                        # X and Y components of vectors weighted by polarization fraction (pfrac, 0)    
                        lenX, lenY = polarization_frac[::SPACING,::SPACING] * np.ones_like(X), np.zeros_like(Y)
                    else:
                        # X and Y components of vectors ( all (1,0) )    
                        lenX, lenY = np.ones_like(X), np.zeros_like(Y)
    
                    x_ticks = np.linspace(0,BOXSIZE_PC,6)
                    y_ticks = x_ticks
                    x_ticks_pos = np.linspace(0,res,6)
                    y_ticks_pos = x_ticks_pos     
                
                    # Plot density with stokes angle vector field overlay 
                    # where vector lengths are polarization fraction
                    print "plotting density_%s" %axis
                    plt.clf()
                    plt.imshow(N, origin='lower', interpolation='nearest', cmap=plt.cm.Blues) #cmap="prism_r"
                    plt.clim(.18,4)
                    plt.colorbar()
                    plt.xticks(x_ticks_pos, x_ticks)
                    plt.yticks(y_ticks_pos, y_ticks)
                    plt.xlabel("size [pc]")                
                    plt.quiver(X, Y, lenX, lenY, angles=theta_stokes[::SPACING,::SPACING],
                               headlength=0.01, headwidth=1, color='red')
         
                    if weight_pfrac:
                        plt.title("Density %s with Polarization Field Pfrac \nWeighted n0 = %d, p = %d frame %d %s" %(axis, param[0], param[1], frame, beta))
                        plt.savefig("%s/density_with_P_field_%s_%d_Pfrac_n0-%04d_p-%d_%s_frame%04d" %(plot_loc, axis, res, param[0], param[1], beta, frame))
                    else:
                        plt.title("Density %s with Polarization Field \nn0 = %d, p = %d frame %d %s" %(axis, param[0], param[1], frame, beta))
                        plt.savefig("%s/density_with_P_field_%s_%d_n0-%d_p-%d_%s_frame%04d" %(plot_loc, axis, res, param[0], param[1], beta, frame)) 

                    print "plotting density with B field%s" %axis
                    plt.clf()
                    plt.imshow(N, origin='lower', interpolation='nearest', cmap=plt.cm.Blues) #cmap="prism_r"
                    plt.clim(.18,4)
                    plt.colorbar()
                    plt.xticks(x_ticks_pos, x_ticks)
                    plt.yticks(y_ticks_pos, y_ticks)
                    plt.xlabel("size [pc]")                
                    plt.quiver(X, Y, lenX, lenY, angles=B_theta[::SPACING,::SPACING],
                               headlength=0.01, headwidth=1, color='red')
                    plt.title("Density %s with Magnetic Field Without Removing Degeneracy frame %d %s" %(axis, frame, beta))
                    plt.savefig("%s/density_with_B_field_without_removing_degeneracy_%s_%d_%s_frame%04d" %(plot_loc, axis, res, beta, frame))                     
                    
                    '''
                    # Plot density with polarization vectors zoomed in on individual leaves
                    if res == 512:
                        SPACING = 2   # pixels between each vector
                    elif res == 8192:
                        SPACING = 4
        
                    Y, X = np.mgrid[0:2*r_lg_max:SPACING, 0:2*r_lg_max:SPACING]
                    #lenX, lenY = polarization_frac[::SPACING,::SPACING] * np.ones_like(X), np.zeros_like(Y)
                    #lenX, lenY = np.ones_like(X), np.zeros_like(Y)   
                    
                    x_ticks = np.abs(np.linspace(-r_lg_pc_max,r_lg_pc_max,6))
                    y_ticks = x_ticks
                    x_ticks_pos = np.linspace(0,2*r_lg_max,6)
                    y_ticks_pos = x_ticks_pos
    
                    for leaf in dg.leaves:
                        print "\nleaf", leaf.idx                        
                        center = leaf.get_peak()[0]
                        ymin = center[0] - r_lg_max
                        ymax = center[0] + r_lg_max
                        xmin = center[1] - r_lg_max
                        xmax = center[1] + r_lg_max  
                        print "ymin,ymax,xmin,xmax = ", ymin,ymax,xmin,xmax
    
                        # Skip leaves whose large disk is beyond domain boundaries
                        if ymin < 0 or xmin < 0 or ymax > res or xmax > res:
                            print 'skip leaf', leaf.idx
                            continue    
                    
                        lenX, lenY = polarization_frac[ymin:ymax:SPACING,xmin:xmax:SPACING] * np.ones_like(X), np.zeros_like(Y)
                                   
                        # make circles around small and large disk                       
                        plt.clf()  
                        plt.close()                      
                        fig, ax = plt.subplots()
                        disk_sm = plt.Circle((r_lg_max,r_lg_max),r_sm,color='green', fill=False)
                        disk_lg = plt.Circle((r_lg_max,r_lg_max),r_lg_max,color='green', fill=False)
                        plt.imshow(N[ymin:ymax, xmin:xmax], origin='lower', interpolation='nearest', cmap=plt.cm.Blues) #cmap="prism_r"
                        ax.autoscale(False)
                        ax.add_patch(disk_sm)
                        ax.add_patch(disk_lg)
                        plt.xticks(x_ticks_pos, x_ticks)
                        plt.yticks(y_ticks_pos, y_ticks)
                        plt.xlabel("radius [pc]")
                        plt.clim(.18,4) 
                        plt.colorbar()
                        plt.quiver(X, Y, lenX, lenY, angles=theta_stokes[ymin:ymax:SPACING,xmin:xmax:SPACING], 
                                   headlength=0.01, headwidth=1, color='red')
                        plt.title("Leaf idx: %d Polarization Map %s %s n0 = %d, p = %d frame %d" %(leaf.idx, axis, beta, param[0], param[1], frame))
                        plt.savefig("%s/leaf%03d_densityImage_%s_n0-%04d_p-%d_%s_frame%04d" %(plot_loc, leaf.idx, axis, param[0], param[1], beta, frame))
                        '''              

def quantify_smoothness_theta_vs_r():
    """
    Quantify the smoothness of the theta(r) plots by calculating
    the variance of theta(r) for each leaf and then take the mean for
    each field strength and axis. Quantify the smoothness another way
    by summing the magnitude of the gradient of theta(r) for each leaf
    and then taking the mean for each field strength and axis.
    """
    file = open("/scratch2/cdb09f/enzotest/%s_res_plots/%s_res_theta_vs_r_smoothness.txt"%(res,res), 'w') 
    for frame in frames:
        file.write("\n\n################### Frame %d ###################\n"%frame)
        file.write("\nMean variance of Theta(r)\n")
        file.write("------------------------- \n\n")
        file.write("         {0:20}{1}".format("Polarization Angle","Mangnetic Field Angle\n"))
        param = [1,1]
        for axis in axes:
            file.write(axis+" axis:\n")
            for sim in sims:
    
                beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]   
    
                at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                plot_loc_leaves = at.get_plot_location_leaves()
                data_loc = at.get_data_location()
    
                P_leaves_delta_theta = bload("%s/P_leaves_delta_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                B_leaves_delta_theta = bload("%s/B_leaves_delta_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
    
                P_variances = []
                B_variances = []
                for angles in P_leaves_delta_theta:
                    P_variances.append(np.var(angles))
                for angles in B_leaves_delta_theta:
                    B_variances.append(np.var(angles))  
    
                avg_P_variance = round(np.mean(P_variances), 2)
                avg_B_variance = round(np.mean(B_variances), 2)
    
                file.write("{0:<4}{1:>17.2f}{2:>20.2f}".format(beta+":", avg_P_variance, avg_B_variance))
                file.write("\n")

        file.write("\n\nMean sum of |dTheta/dr|\n")
        file.write("-----------------------\n\n")
        file.write("         {0:20}{1}".format("Polarization Angle","Mangnetic Field Angle\n"))        
        for axis in axes:
            file.write(axis+" axis:\n")
            for sim in sims:
    
                beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]   
    
                at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                plot_loc_leaves = at.get_plot_location_leaves()
                data_loc = at.get_data_location()
    
                P_leaves_delta_theta = bload("%s/P_leaves_delta_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                B_leaves_delta_theta = bload("%s/B_leaves_delta_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                
                P_grad_sum = []
                B_grad_sum = []
                for angles in P_leaves_delta_theta:
                    P_grad_sum.append(np.sum(np.abs(np.gradient(angles))))
                for angles in B_leaves_delta_theta:
                    B_grad_sum.append(np.sum(np.abs(np.gradient(angles))))
    
                avg_P_grad_sum= round(np.mean(P_grad_sum), 2)
                avg_B_grad_sum= round(np.mean(B_grad_sum), 2)
    
                file.write("{0:<4}{1:>17.2f}{2:>20.2f}".format(beta+":", avg_P_grad_sum, avg_B_grad_sum))
                file.write("\n")
    file.close()            

def make_leaf_contours_plots():
    """
    Outline contours of leaves in density images
    """ 
    for frame in frames:
        # param has the form [n0,p]
        for param in [[1,1]]: 
            #for axis in ['x']:        
            for axis in axes:  
                #for sim in ['High_512']:
                for sim in sims:
                    beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]                
                    at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                    plot_loc = at.get_plot_location_leaf_contours()
                    data_loc = at.get_data_location()
    
                    N = at.get_fits_array('density',axis)
                    d = Dendrogram.compute(N, min_value=2.0, min_delta=1.0, min_npix=5, verbose=True)
                    p = d.plotter()    
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 1, 1)
                    ax.imshow(N, origin='lower', interpolation='nearest', cmap=plt.cm.Blues, vmax=4.0, vmin=0.18)
                    plt.title("Leaves %s %s frame %d"%(axis, beta, frame))
                    for leaf in d.leaves:
                        p.plot_contour(ax, structure=leaf.idx, lw=0.5, colors='red')
                    fig.savefig("%s/leaf_contours_%s_%s_frame%04d"%(plot_loc, axis, beta, frame))

def make_leaves_html_tables():
    """
    make html tables with leaf image on lthe eft and leaf's theta(r) plot on the right.
    one table for each axis, field strength, and frame.
    """
    param = [1,1]
    for axis in axes:
        for sim in sims:
            beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]
            for frame in frames:

                at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                data_loc = at.get_data_location()
                leafidx = bload("%s/leafidx_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                file_name = at.get_plot_location_leaves() + '/leaves_table_%s_%s_%s_frame%04d.html'%(res,axis,beta,frame)
                
                file = open(file_name, 'w')
                file.write("<html>\n")
                file.write("<table>\n")

                for idx in leafidx:
                    image = "leaf%03d_densityImage_%s_n0-0001_p-1_%s_frame%04d.png"%(idx,axis,beta,frame)
                    graph = "leaf%03d_P_delta_theta_vs_r_disk_%s_n0-0001_p-1_%s_frame%04d.png"%(idx,axis,beta,frame)
                    
                    file.write("  <tr>\n")
                    file.write("    <td><img src=%s></td>\n"%image)
                    file.write("    <td><img src=%s></td>\n"%graph)
                    file.write("  </tr>\n")

                file.write("</table>\n")
                file.write("</html>\n")
                file.close()

def make_mass_vs_smoothness_plots():
    """
    plot mass of leaf vs variance of delta_theta(r), mass in small disk vs variance of delta_theta(r)
    make similar plots but with sum of |d_delta_theta(r)/dr| instead of variance
    make plots of mass in small or large disk vs delta_theta
    """
    param = [1,1]
    for frame in [50]:
    #for frame in frames:
        for axis in axes:
            for sim in sims:
                beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]
                at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                data_loc = at.get_data_location()
                plot_loc = at.get_plot_location_mass_vs_smoothness()

                leafidx = bload("%s/leafidx_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                P_leaves_delta_theta = bload("%s/P_leaves_delta_theta_whole_box_%s_n0-%04d_p-%d_frame%04d_%s.data" %(data_loc, axis, param[0], param[1], frame, beta))
                B_leaves_delta_theta = bload("%s/B_leaves_delta_theta_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                mass_in_sm_disk = bload("%s/mass_in_sm_disk_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))
                mass_leaves = bload("%s/mass_leaves_whole_box_%s_frame%04d_%s.data" %(data_loc, axis, frame, beta))           
 
                #dg = at.get_dendrogram('density', axis)

                P_delta_theta = np.array([item[list(r_disk).index(r_lg_pc_max)] for item in P_leaves_delta_theta])
                mass_in_lg_disk = np.array([item[list(r_disk).index(r_lg_pc_max)] for item in mass_leaves])
                
                mass = []
                P_variances = []
                B_variances = []
                P_grad_sum = []
                B_grad_sum = []               

                #for idx in leafidx:
                #    mass.append(np.sum(dg[idx].values()))
                #    if np.sum(dg[idx].values()) > 2000:
                #        print "leaf %d: mass = %d, n_pix = %d, frame %d, %s, %s"%(idx, round(np.sum(dg[idx].values())), dg[idx].get_npix(), frame, axis, beta)
                
                P_variances = [np.var(angles) for angles in P_leaves_delta_theta]
                B_variances = [np.var(angles) for angles in B_leaves_delta_theta]
                P_grad_sum = [np.sum(np.abs(np.gradient(angles))) for angles in P_leaves_delta_theta]
                B_grad_sum = [np.sum(np.abs(np.gradient(angles))) for angles in B_leaves_delta_theta]

                mass = np.array(mass)
                P_variances = np.array(P_variances)
                B_variances = np.array(B_variances)
                P_grad_sum = np.array(P_grad_sum)
                B_grad_sum = np.array(B_grad_sum)

                ## Scatter plots smoothness vs leaf mass
                #print" len(mass), len(P_variances)",  len(mass), len(P_variances)
                #plt.clf()
                #plt.title('Leaf Mass vs Polarization Angle Variance\nframe %d, %s, %s'%(frame, axis, beta))
                #plt.xlabel("Leaf Mass")
                #plt.ylabel(r'$\sigma_{\theta}^2$')              
                #plt.xlim(0,mass.max())
                #plt.ylim(0,P_variances.max())
                #plt.scatter(mass, P_variances)
                #plt.savefig("%s/P_theta_variance_vs_leaf_mass_frame%04d_n0-%04d_p-%d_%s_%s_" %(plot_loc, frame, param[0], param[1], axis, beta))                
                #
                #"""
                #plt.clf()
                #plt.title('Leaf Mass vs Magnetic Field Angle Variance\nframe %d, %s, %s'%(frame, axis, beta))
                #plt.xlabel("Leaf Mass")
                #plt.ylabel(r'$\sigma_{\theta}^2$')              
                #plt.xlim(0,mass.max())
                #plt.ylim(0,B_variances.max())
                #plt.scatter(mass, B_variances, color='black')
                #plt.savefig("%s/B_theta_variance_vs_leaf_mass_frame%04d_%s_%s_" %(plot_loc, frame, axis, beta)) 
                #"""
#
                #plt.clf()
                #print" len(mass), len(P_grad_sum)",  len(mass), len(P_grad_sum)
                #plt.title('Leaf Mass vs Polarization Angle '+r'$\Sigma d\theta/dr$'+'\nframe %d, %s, %s'%(frame, axis, beta))
                #plt.xlabel("Leaf Mass")
                #plt.ylabel(r'$\Sigma d\theta/dr$')              
                #plt.xlim(0,mass.max())
                #plt.ylim(0,P_grad_sum.max())
                #plt.scatter(mass, P_grad_sum)
                #plt.savefig("%s/P_sum-grad-theta_vs_leaf_mass_frame%04d_n0-%04d_p-%d_%s_%s_" %(plot_loc, frame, param[0], param[1], axis, beta))
                #
                #"""
                #plt.clf()
                #plt.title('Leaf Mass vs Magnetic Field Angle '+r'$\Sigma d\theta/dr$'+'\nframe %d, %s, %s'%(frame, axis, beta))
                #plt.xlabel("Leaf Mass")
                #plt.ylabel(r'$\Sigma d\theta/dr$')              
                #plt.xlim(0,mass.max())
                #plt.ylim(0,B_grad_sum.max())
                #plt.scatter(mass, B_grad_sum, color='black')
                #plt.savefig("%s/B_sum-grad-theta_vs_leaf_mass_frame%04d_%s_%s_" %(plot_loc, frame, axis, beta))
                #"""

                # Scatter plots smoothness vs mass in disks
                plt.clf()
                plt.title('Mass in Small Disk vs Polarization Angle Variance\nframe %d, %s, %s'%(frame, axis, beta))
                plt.xlabel("Mass in Small Disk")
                plt.ylabel(r'$\sigma_{\theta}^2$')              
                #plt.xlim(mass_in_sm_disk.min(),mass_in_sm_disk.max())
                plt.xlim(0,mass_in_sm_disk.max())
                plt.ylim(0,P_variances.max())
                plt.scatter(mass_in_sm_disk, P_variances)
                plt.savefig("%s/P_theta_variance_vs_small_disk_mass_frame%04d_n0-%04d_p-%d_%s_%s_" %(plot_loc, frame, param[0], param[1], axis, beta))

                plt.clf()
                plt.title('Mass in Small Disk vs Polarization Angle '+r'$\Sigma |d\theta/dr|$'+'\nframe %d, %s, %s'%(frame, axis, beta))
                plt.xlabel("Mass in Small Disk ")
                plt.ylabel(r'$\Sigma d\theta/dr$')              
                #plt.xlim(mass_in_sm_disk.min(),mass_in_sm_disk.max())
                plt.xlim(0,mass_in_sm_disk.max())
                plt.ylim(0,P_grad_sum.max())
                plt.scatter(mass_in_sm_disk, P_grad_sum)
                plt.savefig("%s/P_sum-grad-theta_vs_small_disk_mass_frame%04d_n0-%04d_p-%d_%s_%s_" %(plot_loc, frame, param[0], param[1], axis, beta))                               

                plt.clf()
                print" len(mass_in_sm_disk), len(P_delta_theta)",  len(mass_in_sm_disk), len(P_delta_theta)
                print "max(mass_in_sm_disk), max(P_delta_theta)", max(mass_in_sm_disk), max(P_delta_theta)
                plt.title('Mass in Small Disk vs '+r'$\Delta\theta$'+' Polarization\nframe %d, %s, %s'%(frame, axis, beta))
                plt.xlabel("Mass in Small Disk ")
                plt.ylabel(r'$\Delta\theta$', size=15)              
                #plt.xlim(mass_in_sm_disk.min(),mass_in_sm_disk.max())
                plt.xlim(0,mass_in_sm_disk.max())
                plt.ylim(0,P_delta_theta.max())
                plt.scatter(mass_in_sm_disk, P_delta_theta)
                plt.savefig("%s/P_delta_theta_vs_small_disk_mass_frame%04d_n0-%04d_p-%d_%s_%s_" %(plot_loc, frame, param[0], param[1], axis, beta)) 

                plt.clf()
                print" len(mass_in_lg_disk), len(P_delta_theta)",  len(mass_in_lg_disk), len(P_delta_theta)
                print "max(mass_in_lg_disk), max(P_delta_theta)", max(mass_in_lg_disk), max(P_delta_theta)
                plt.title('Mass in Small Disk vs '+r'$\Delta\theta$'+' Polarization\nframe %d, %s, %s'%(frame, axis, beta))
                plt.xlabel("Mass in Small Disk ")
                plt.ylabel(r'$\Delta\theta$', size=15)              
                #plt.xlim(mass_in_lg_disk.min(),mass_in_lg_disk.max())
                plt.xlim(0,mass_in_lg_disk.max())
                plt.ylim(0,P_delta_theta.max())
                plt.scatter(mass_in_lg_disk, P_delta_theta)
                plt.savefig("%s/P_delta_theta_vs_large_disk_mass_frame%04d_n0-%04d_p-%d_%s_%s_" %(plot_loc, frame, param[0], param[1], axis, beta))                 
                
def mean_B_theta_whole_box():
    param=[1,1]
    for frame in frames:
        print "\n"+str(frame)
        for axis in axes:
            print "\n"
            for sim in sims:
                beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]
                field_horizontal = {'x':'By','y':'Bz','z':'Bx'}[axis]
                field_vertical   = {'x':'Bz','y':'Bx','z':'By'}[axis] 
                at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
                aB_hor = at.get_fits_array(field_horizontal, axis)
                aB_vert= at.get_fits_array(field_vertical, axis)  
                Btheta_boxAvg = np.mean(np.arctan2(aB_vert,aB_hor))*180/np.pi

                print "Btheta_boxAvg", axis, beta, np.round(Btheta_boxAvg, 2)

def weird_corners_theta():
    param = [1,1]
    frame = 50
    axis = 'y'
    for sim in sims:
        beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]
        field_horizontal = {'x':'By','y':'Bz','z':'Bx'}[axis]
        field_vertical   = {'x':'Bz','y':'Bx','z':'By'}[axis]         
        at = access_thing(sim,res,frame,param[0],param[1], machine=machine)

        # Get frbs
        Q = at.get_fits_array('Q%s_n0-%04d_p-%d'%(axis,param[0],param[1]), axis)
        U = at.get_fits_array('U%s_n0-%04d_p-%d'%(axis,param[0],param[1]), axis)
        B_hor = at.get_fits_array(field_horizontal, axis)
        B_vert= at.get_fits_array(field_vertical, axis)

        P_theta = 0.5*np.arctan2(U,Q)*180/np.pi 
        B_theta = np.arctan2(B_vert, B_hor)*180/np.pi 
        #B_theta[B_theta > 90] = B_theta[B_theta > 90] - 180
        #B_theta[B_theta < -90] = 180 + B_theta[B_theta < -90]
        
        two_pc_pix = 2.0 * PC_TO_PIXEL

        plt.clf()
        fig, ax = plt.subplots()
        plt.title("Polarization Angle Image Y %s" %beta)
        plt.imshow(P_theta, origin='lower', interpolation='nearest', cmap=plt.cm.Greys)
        plt.colorbar() 
        disk_lg = plt.Circle((res/2,res/2),two_pc_pix,color='red', fill=False)
        ax.autoscale(False)
        ax.add_patch(disk_lg)
        plt.savefig('theta_images/P_theta_image_y_disk2pc_%s'%beta) 

        plt.clf()
        fig, ax = plt.subplots()
        plt.title("B Field Angle Image Without Removing Degeneracy %s" %beta)
        plt.imshow(B_theta, origin='lower', interpolation='nearest', cmap=plt.cm.Greys) 
        plt.colorbar()        
        disk_lg = plt.Circle((res/2,res/2),two_pc_pix,color='red', fill=False)
        ax.autoscale(False)
        ax.add_patch(disk_lg)
        plt.savefig('theta_images/B_theta_image_without_removing_degeneracy_y_disk2pc_%s'%beta) 

def make_test_data():
    B0 = 1
    Bx = np.ones([512.,512.,512.])
    By = Bx
    Bz = -Bx
    fake1 = {'Bx':Bx, 'By':By, 'Bz':Bz}
    Q_map_fake = _Qlocal()

def CMB_E_and_B_modes():
    """
    Decompose polarazation map (Stokes Q and U) into E and B modes. Calculate the ratio of E to B.
    E is the curl free component of the vector field
    B is the div free component of the vector field
    """
    res = 512
    frame = 50
    param = [1,1]

    res = 512
    N = np.array([res,res], dtype=np.int32)
    xsize = 1 * np.pi / 180
    size2d = np.array([xsize,xsize])
    Delta = size2d/N
    print "Delta", Delta
    Deltal = cmbtools.Delta2l(Delta, N)
    print "Deltal", Deltal

    for axis in axes:
        for sim in sims:
            beta = {'High_512':'B02', 'Mid_512':'B2', 'Low_512':'B20'}[sim]       
            at = access_thing(sim,res,frame,param[0],param[1], machine=machine)
            plotloc = at.get_plot_location_QU2EB()
    
            # Get frbs
            Q = np.array(at.get_fits_array('Q%s_n0-%04d_p-%d'%(axis,param[0],param[1]), axis))
            U = np.array(at.get_fits_array('U%s_n0-%04d_p-%d'%(axis,param[0],param[1]), axis))

            print "Q\n", Q 
            print "U\n", U
            print "Q.shape", Q.shape
            print "U.shape", U.shape

            # 2d fourier transforms of Q and U
            #Qharm = cmbtools.map2harm(Q, Delta)
            #Uharm = cmbtools.map2harm(U, Delta)
            Qharm = np.fft.fft2(Q)
            Uharm = np.fft.fft2(U)

            print "Qharm\n", Qharm
            print "Uharm\n", Uharm
            
            # Decompose stokes into E and B modes
            Eharm, Bharm = cmbtools.QU2EB(Qharm, Uharm, Deltal)

            print "Eharm\n", Eharm
            print "Bharm\n", Bharm
            
            # 2d inverse fourier transforms of E and B to bring them back into real space
            #E = cmbtools.harm2map(Eharm, Delta)
            #B = cmbtools.harm2map(Bharm, Delta)
            E = np.fft.ifft2(Eharm)
            B = np.fft.ifft2(Bharm)

            print "E\n", E
            print "B\n", B

            # Compute coeffictients for the harmonics
            deltalbins=2.5
            lmax = 100
            lbins = np.arange(0,lmax,deltalbins)

            ClE = cmbtools.harm2cl(Eharm, Deltal, lbins)
            ClB = cmbtools.harm2cl(Bharm, Deltal, lbins)
            print "ClE/ClB", ClE/ClB

            plt.clf()
            plt.imshow(Q,interpolation='nearest')
            plt.title('Q')
            plt.colorbar()
            plt.savefig('p49_Q_%s_%s'%(axis, beta))
            
            plt.clf()
            plt.imshow(U,interpolation='nearest')
            plt.title('U')
            plt.colorbar()
            plt.savefig('p49_U_%s_%s'%(axis, beta))

            plt.clf()
            plt.imshow(E,interpolation='nearest')
            plt.title('E')
            plt.colorbar()
            plt.savefig('p49_E_%s_%s'%(axis, beta))
            
            plt.clf()
            plt.imshow(B,interpolation='nearest')
            plt.title('B')
            plt.colorbar()
            plt.savefig('p49_B_%s_%s'%(axis, beta))
            
         

if 0:
    make_frbs()
if 0:
    analyze_dense_regions()

if 0:
    print "making P_frac plots"
    sys.stdout.flush()
    make_P_frac_plots()
if 0:
    print "making field strength plots"
    sys.stdout.flush()
    make_field_strength_plots()
if 0:
    print "making viewing angle plots"
    sys.stdout.flush()
    make_viewing_angle_plots()
if 0:
    print "making threshold density plots"
    sys.stdout.flush()
    make_threshold_density_plots()
if 0:
    print "making time dependence plots"
    sys.stdout.flush()
    make_time_dependence_plots()     
if 0:
    print "making Pol vs B theta plots"
    sys.stdout.flush()
    make_P_delta_theta_vs_B_delta_theta_plots()
if 0:
    print "making density with B field plots"
    sys.stdout.flush()
    make_density_with_Bfield_plots()
if 0:
    print "making theta_vs_r plots"
    sys.stdout.flush()
    make_theta_vs_r_plots()       
if 0:
    print "calculating theta(r) smoothness"
    sys.stdout.flush()
    quantify_smoothness_theta_vs_r()
if 0:
    print "make leaf contours plots"    
    sys.stdout.flush()
    make_leaf_contours_plots()
if 0:
    print "make leaves html tables"
    sys.stdout.flush()
    make_leaves_html_tables()
if 0:
    print "make mass vs smoothness plots"
    sys.stdout.flush()
    make_mass_vs_smoothness_plots()
if 0:
    weird_corners_theta()
if 1:
    CMB_E_and_B_modes()


endtime = time.time() - start
print "\n\nFINISHED\n\n"
print "elapsed time = %d min %d sec\n" %(endtime/60, endtime%60)

    
                       
