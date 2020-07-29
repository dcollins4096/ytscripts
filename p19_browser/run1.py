
from go import *
import product
reload(product)
import make_page
reload(make_page)
basedir='/Users/dcollins/RESEARCH3/Paper19_47_overlap/0000_main_plots/'
glob1 = "%s/density_radius/density_radius_c????.png"%basedir
reg1  = re.compile(r"%s/density_radius/density_radius_c(\d\d\d\d).png"%basedir)
p1 = product.product('rho_r',regexp=reg1,myglob=glob1,parameters=['core_id'],style='single')
p1.get_frames()

#glob1 = "%s/core_proj_follow/density_radius_c????.png"%basedir
#glob2  = r"%s/core_proj_follow/follow_c????_n????_centered_Projection_x_density.png"%(basedir)
#reg2   = re.compile(r"%s/core_proj_follow/follow_c(\d\d\d\d)_n(\d\d\d\d)_centered_Projection_x_density.png"%(basedir))
#p2 = product.product('core_proj_follow',regexp=reg2,myglob=glob2,parameters=['core_id','frame'],style='frames')
##p2.check_glob()
#p2.get_frames()

cl=make_page.make_page([p1,p2])
