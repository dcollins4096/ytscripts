
from go import *

loc="/scratch/00369/tg456484/g47_xsede_2022/Turbulence"
ol = loc+"/meta_%d/OutputLog"

mach_list = [10,20,40]

dicts={}

for mach in mach_list:
    fptr = open(ol%mach)
    lines=fptr.readlines()
    fptr.close()
    names = []
    cycle = []
    time = []
    wall = []
    for line in lines:
        things = line.split()
        names.append(things[2])
        cycle.append(float(things[3]))
        time.append(float(things[4]))
        wall.append(float(things[5]))

    stuff = {'names':names,'cycle':nar(cycle),'time':nar(time),'wall':nar(wall)}
    start_cycle = -10
    stop_cycle = -1
    DN = cycle[stop_cycle]-cycle[start_cycle]
    DT = time[stop_cycle]-time[start_cycle]
    dt_sim = DT/DN
    DWall = wall[stop_cycle]-wall[start_cycle]
    dt_wall = DWall/DN
    Nzones = 256**3
    Nnodes=1

    seconds_update = dt_wall
    seconds_zone_up = dt_wall/Nzones
    hour_zone_up = seconds_zone_up/3600
    SU_zone_up = hour_zone_up*Nnodes
    mach_number= float(mach)/10
    fiducial_mach = 1

    print(DT, dt_sim, dt_wall, SU_zone_up, dt_sim*(1+mach_number)/(1+fiducial_mach))

    dicts[mach]=stuff


