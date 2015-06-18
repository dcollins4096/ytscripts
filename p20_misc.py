import pyfits
from yt.analysis_modules.level_sets.api import * #for clumps
import pyximport; pyximport.install()
import particle_ops
import random
import clump_particles
from yt.utilities.data_point_utilities import FindBindingEnergy
from yt.utilities.physical_constants import \
            gravitational_constant_cgs as G
reload(clump_particles)
ef('particle_selector.py')
def annotate_box(pw,L,R,plot_args):
    pw.annotate_image_line(L,[R[0],L[1],1], plot_args=plot_args)
    pw.annotate_image_line(L,[L[0],R[1],1], plot_args=plot_args)
    pw.annotate_image_line([L[0],R[1],0],R, plot_args=plot_args)
    pw.annotate_image_line([R[0],L[1],0],R, plot_args=plot_args)

def write_fits_2d(dirname, frame, simname, axis, field,resolution, dir_2d):
    ds = yt.load('%s/DD%04d/data%04d'%(dirname,frame,frame))
    proj = ds.proj(field,axis)
    frb =  proj.to_frb(1.0,resolution)
    hdu = pyfits.PrimaryHDU(frb[field])
    set_2d = '%s/%s_n%04d_r%04d_%s_%s.fits'%(dir_2d, simname, frame, resolution,'xyz'[axis], field)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(set_2d)
    hdulist.close()
    
def write_fits_3d(dirname, frame, simname, axis, field,resolution,level, dir_3d):
    """ *resolution* is a single element, though covering grid wants a tripple. Assumes cubic."""
    ds = yt.load('%s/DD%04d/data%04d'%(dirname,frame,frame))
    cg = ds.covering_grid(level,[0.0]*3, [resolution]*3)
    hdu = pyfits.PrimaryHDU(cg[field])
    set_3d = '%s/%s_n%04d_r%04d_%s_%s.fits'%(dir_3d, simname, frame, resolution,'3d', field)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(set_3d)
    hdulist.close()
    print "wrote", set_3d

def plot_disperse(set_density, set_disperse, figurename, dirname1, dirname2=None, filament_list=None):
    if dirname2 is None:
        dirname2 = dirname1
    base_figname, base_format = tuple(figurename.split('.'))
    density_fits = pyfits.open("%s/%s"%(dirname1,set_density))
    density_full = density_fits[0].data
    density_fits.close
    disp_fits = pyfits.open("%s/%s"%(dirname2,set_disperse))
    disp_full = disp_fits[0].data
    disp_fits.close
    rank = len(disp_full.shape)
    if filament_list is None:
        fils = np.unique(disp_full)
        all_or_one = 'all'
    else:
        fils = filament_list
        all_or_one = filament_list[0]
    nfils = len(fils)
    rmap = rainbow_map(nfils)
    right_x = disp_full.shape[0]
    right_y = disp_full.shape[1]
    dx=1
    y, x = np.mgrid[0.5*dx:right_x-0.5*dx:right_x*1j, 0.5*dx:right_y-0.5*dx:right_y*1j]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if rank == 2:
        z_coords = [0]
    else:
        z_coords = range(disp_full.shape[2])
        z_coords = [0,1,2]
    for z in z_coords:
        if rank == 2: 
            density = density_full
            disp = disp_full
            this_fig_name = figurename
            this_fig_name = "%s_fil_%s.%s"%(base_figname,all_or_one, base_format)
            density_fig_name = "%s_den.%s"%(base_figname, base_format)

        else:
            density = density_full[:,:,z]
            disp = disp_full[:,:,z]
            this_fig_name = "%s_z%04d_fil_%s.%s"%(base_figname,z,all_or_one, base_format)
            density_fig_name = "%s_z%04d_den.%s"%(base_figname, z, base_format)
        ax.cla()
        ax.imshow( np.log10( density) ,cmap='gray', origin='lower')
        fig.savefig(density_fig_name)
        print "saved", density_fig_name

        if filament_list is None:
            these_fils = np.unique(disp)[1:]
        else:
            these_fils = filament_list
        for nfil in these_fils:
            this_fil = disp==nfil
            if this_fil.size > 0:
                #print " ", nfil,
                x_fil = x[this_fil]
                y_fil = y[this_fil]
                ax.scatter(x_fil,y_fil, c=rmap(nfil))
                avg_x = x_fil.sum()/x_fil.size
                avg_y = y_fil.sum()/y_fil.size
                ax.text(avg_x,avg_y,nfil,color=rmap(nfil), fontsize=8)

        fig.savefig(this_fig_name)
        print "saved", this_fig_name
    plt.close(fig)
    return density_full, disp_full

#def get_3d_from_2d(fil_2d,disp_2d,disp_3d)
def bd(value, dx):
    """bump down"""
    return np.floor( value/dx )*dx
def bu(value, dx):
    """bump up"""
    return np.floor( value/dx )*dx
def get_region(dirname, frame, simname, axis, fil_id, disp_2d,prefix,x_range=[0.0,1.0]):
    prefix2 = prefix+"_fil%04d"%fil_id
    right_x = disp_2d.shape[0]
    right_y = disp_2d.shape[1]
    
    dx=1.
    delta_x = dx/right_x
    y, x = np.mgrid[0.5*dx:right_x-0.5*dx:right_x*1j, 0.5*dx:right_y-0.5*dx:right_y*1j]/right_x
    this_fil = disp_2d==fil_id
    x_fil = x[this_fil]
    y_fil = y[this_fil]
    fil_left = bd(x_fil.min() - 3*delta_x, delta_x), bd(y_fil.min()-3*delta_x, delta_x)
    reg_left   =  [ x_range[0], fil_left[0], fil_left[1] ]
    fil_right  =  [ bu(x_fil.max() + 3*delta_x, delta_x), bu(y_fil.max()+3*delta_x, delta_x)]
    reg_right  =  [ x_range[1], fil_right[0], fil_right[1]]
    reg_left=nar(reg_left); reg_right=nar(reg_right)

    print reg_left
    print reg_right
    #reg_left=np.zeros(3)
    #reg_right=np.ones(3)
    reg_cen = 0.5*(reg_left+reg_right)
    ds = yt.load('%s/DD%04d/data%04d'%(dirname,frame,frame))
    reg = ds.region(reg_cen, reg_left, reg_right)
    return reg_cen, reg_left, reg_right, reg
def transverse(dirname, frame, simname, axis, fil_id, disp_2d,prefix,x_range=[0.0,1.0]):
    ds = yt.load('%s/DD%04d/data%04d'%(dirname,frame,frame))
    reg_cen, reg_left, reg_right, reg = get_region(dirname, frame, simname, axis, fil_id, disp_2d,prefix,x_range)
    reg = ds.region(reg_cen, reg_left, reg_right)
    proj = ds.proj(field,0)
    pw=proj.to_pw()
    pw.set_cmap('density','gray')
    annotate_box(pw, fil_left, fil_right, plot_args={'color':'r'})
    print pw.save(prefix2+"_full")
    for ax in [0,1,2]:
        proj = ds.proj(field,ax, data_source = reg)
        pw=proj.to_pw() 
        pw.set_cmap('density','gray')
        #pw.set_zlim('density',zlim[0],zlim[1])
        print pw.save(prefix2+"sub")
    return reg

def preimage(dirname, target_frame, axis, source_region, prefix,center=[0.5]*3):
    field = 'density'
    ds = yt.load('%s/DD%04d/data%04d'%(dirname,target_frame,target_frame))
    proj = ds.proj(field,axis,center=center)
    pw = proj.to_pw(center=center)
    pw.set_cmap('density','gray')
    source_indices = source_region['particle_index']
    pw.annotate_dave_particles(1.0, indices = source_indices, col='y')
    source_indices = source_region['particle_index']
    outname = "%s_nt%04d"%(prefix, target_frame)
    pw.save(outname)
    print "saved", outname


disp_name_all = 'u05_n0125_r0128_x_density.fits.up.NDskl.fits'
disp_name_cut1 = 'u05_n0125_r0128_x_density.fits_c1.up.NDskl.fits'
disp_name_3d_cut1 = 'u05_n0125_r0128_3d_density.fits_c1.up.NDskl.fits'

dirname = '/scratch1/dcollins/Paper19/u05-r4-l4-128'
frame = 125
#frame = 20
simname = 'u05'
axis = 0
field = 'density'
resolution = 128
dir_2d = '/scratch1/dcollins/Paper37_Filaments/2d/'
dir_3d = '/scratch1/dcollins/Paper37_Filaments/3d/'
cut = 1

set_2d = '%s_n%04d_r%04d_%s_%s.fits'%( simname, frame, resolution,'xyz'[axis], field)
fils_2d = '%s_n%04d_r%04d_%s_%s.fits_c%d.up.NDskl.fits'%( simname, frame, resolution,'xyz'[axis], field,cut)
outname_image = '%s_n%04d_r%04d_%s_%s_PROJ.png'%( simname, frame, resolution,'xyz'[axis], field)
#density_2d, disp_2d = plot_disperse(set_2d,fils_2d,outname_image , dir_2d,filament_list=[90,91])
#reg=transverse(dirname, frame, simname, axis, 90, disp_2d, prefix='u05_n200_reg_xcut_2', x_range=[0.075,0.125])
reg_cen, reg_left, reg_right, reg = get_region(dirname, frame, simname, axis, 90, disp_2d, prefix='u05_n200_reg_xcut_2', x_range=[0.075,0.125])
for target_frame in [0]: #,10,20,30,40,50,60,70,80,90,100,110,120,125]:
    preimage(dirname, target_frame, 0, reg, 'preimage_ns0125_fil_0090')
#density_3d, disp_3d = plot_disperse('u05_n0125_r0128_3d_density.fits',disp_name_3d_cut1, 'test_3d.png', dir_3d)
#write_fits_2d(dirname, frame, simname, axis, field, resolution, dir_2d)
#write_fits_3d(dirname, frame, simname, axis, field, resolution, 0, dir_3d)
#disp3d = pyfits.open("%s/%s"%(dir_3d,'u05_n0125_r0128_3d_density.fits_c1.up.NDskl.fits'))[0].data


