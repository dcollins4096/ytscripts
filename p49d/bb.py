from go import *
import astropy.io.fits as pyfits
if 'p42_new_driving_ak' not in sys.path:
    sys.path.append('p42_new_driving_ak')
if '/scratch1/dcollins/Paper49c/p49_eigenmodes' not in sys.path:
    sys.path.append("/scratch1/dcollins/Paper49c/p49_eigenmodes")
import p42_helmholtz as p42h
reload(p42h)
import taxi_subclass
reload(taxi_subclass)

import p49_fields

#bcar = taxi.load('bb_b1p1')
bcar = taxi.load('bb_b1p.032')
bcar.derived_fields['QU'] = p49_fields.add_QU

all_bb=['bb_b1p.01' , 'bb_b1p.032', 'bb_b1p.1'  , 'bb_b1p.32' , 'bb_b1p1'   ,\
'bb_b.1p.01', 'bb_b.1p.032', 'bb_b.1p.1' , 'bb_b.1p.32', 'bb_b.1p1'  ]
flt = taxi.fleet(all_bb)

if 1:
    import turb_quan
    reload(turb_quan)
    for car in flt:
        car.derived_fields['QU'] = p49_fields.add_QU
        qb = turb_quan.quan_box(car=car)
        qb.EBall()
        qb()

"""


if False:
    if 'data3' not in dir():
        import p49_eigen
        reload(p49_eigen)
        bcar.load(frame)
        data3 = p49_eigen.get_cubes_uniform(bcar.data)
    def check_dir_make(dirname): 
        if not os.path.exists(dirname):
            os.mkdir(dirname)

    if 'data3' in dir():

        import p49_plot_tools
        analysis={#'print_wave':False,
                  'plot_fields':1,
                  'k_mag':1,
                  'k_proj':1,
                  'k_func':p49_plot_tools.maxis
                 }

        plot_dir = './p49_bb_plots/'
        this_plot_dir = "%s/%04d"%(plot_dir,frame)
        check_dir_make(plot_dir)
        prefix='SSS'
        check_dir_make(this_plot_dir)
        p49_plot_tools.do_stuff(data3,outdir="%s/%s_%04d_"%(this_plot_dir,prefix,frame),**analysis)

if 'bbu' not in dir() and False:
    bbu = p42h.short_oober()
    bbu.data = data
    bbu.directory = "/scratch1/dcollins/Paper49_EBQU/bburkhart/MHD/256/b1p1"

"""
