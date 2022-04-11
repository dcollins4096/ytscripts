
from go import *

loc="g47/p49_outputs/"
ol = loc+"/%s/OutputLog"
mach_list = [10,20,40]

dicts={}

sims=["1_1", "1_2", "1_half", "2_1", "2_2", "2_half", "3_1", "3_2", "3_half", "half_1", "half_2", "half_half"]
for sim in sims:
    fptr = open(ol%sim)
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
    Nzones = 512**3
    Nnodes=1



    seconds_update = dt_wall
    seconds_zone_up = dt_wall/Nzones
    hour_zone_up = seconds_zone_up/3600
    SU_zone_up = hour_zone_up*Nnodes
    stuff['DN']=DN
    stuff['DT']=DT
    stuff['dt_sim']=dt_sim
    stuff['DWall']=DWall
    stuff['dt_wall']=dt_wall
    stuff['SU_zone_up']=SU_zone_up
    print(SU_zone_up)
    dicts[sim]=stuff

