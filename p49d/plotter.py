import sys
import os
from mpi4py import MPI
sys.path = ['.']+sys.path
sys.path.append('tools_yt_etc')
from go import *

#sys.path = ["%s/repos/p49d_turb_sims/python"%os.environ['HOME']]+sys.path
#from GL import *
import yt
import taxi
reload(taxi)
yt.config.ytcfg['yt','loglevel']="50"  #keeps it quiet
yt.enable_parallelism()


###
from optparse import OptionParser
usage = "usage: %plotter.py [options] "
parser = OptionParser(usage=usage)
parser.add_option('-p', '--projections', dest='proj', action='store_true', 
                  help='do projections', default=False)
parser.add_option('-q', '--queb', dest='queb', action='store_true',
                  help='Do Q,U,E,B plots', default=False)
parser.add_option('-s', '--stats', dest='stats', action='store_true',
                  help='compute first and second order statistics',default=False )
(options, args) = parser.parse_args()
###

this_car = taxi.load( args[0])
print(options,args)

if options.proj:
    print('ok')
    for frame in this_car.return_frames():
        ds = this_car.load(frame)
        for ax in 'xyz': #this_car.axis:
            for field in this_car.fields:
                "data0022_Projection_x_density.png"
                output_name = "%s_%04d_Projection_%s_%s.png"%(this_car.outname,frame,ax,field)
                if os.path.exists(output_name):
                    print("Plot exists, continue %s"%output_name)
                    continue
                plot=yt.ProjectionPlot(ds,ax,field)
                if yt.is_root():
                    plot.save(output_name)

import time
comm = MPI.COMM_WORLD
rank = comm.Get_rank()



if options.stats:
    import turb_quan
    reload(turb_quan)
    qb = turb_quan.quan_box(this_car)
    qb.load()
    tstart = time.time()
    qb()
    tend = time.time()
#   if rank == 0:
#       fptr = open("mpi_test.txt","a")
#       fptr.write("=== allquan nproc %d  frame %d === \n"%(comm.Get_size(), frame))
#       fptr.write("tstart %0.16f\n"%tstart)
#       fptr.write("tstop  %0.16f\n"%tend)
#       fptr.write("tdiff %0.16f\n"%(tend-tstart))
#       fptr.close()


if options.queb:
    import turb_quan
    reload(turb_quan)
    qb = turb_quan.quan_box(this_car)
    for frame in this_car.return_frames():
        tstart = time.time()
        qb.make_frbs(frame)
        tend = time.time()
        if rank == 0:
            fptr = open("mpi_test.txt","a")
            fptr.write("=== frbs nproc %d  frame %d === \n"%(comm.Get_size(), frame))
            fptr.write("tstart %0.16f\n"%tstart)
            fptr.write("tstop  %0.16f\n"%tend)
            fptr.write("tdiff %0.16f\n"%(tend-tstart))
            fptr.close()
    #qb.EBall()



