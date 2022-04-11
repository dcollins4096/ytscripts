from go import *

probably did not use this code.

arr=nar([[0.1 , 5.00   , 46.9 , 11.0 , 41.8 , 10.7 , 13.3 ],
[0.5 , 1.00   , 26.5 , 16.4 , 15.7 , 3.77 , 1.76 ],
[1.0 , 0.500  , 9.62 , 14.9 , 14.2 , 5.20 , 2.20 ],
[2.0 , 0.250  , 14.2 , 14.5 , 13.8 , 1.76 , 0.723],
[4.0 , 0.125  , 4.15 , 14.6 , 13.9 , 3.01 , 1.25 ],
[7.0 , 0.0714 , 2.06 , 14.5 , 13.9 , 3.47 , 1.43 ]])

s1=slice(1,None)
Mach = arr[s1,0]
Tdyn = arr[s1,1]
dt_sim=arr[s1,2]*1e-6
dt_wal=arr[s1,3]
mus_zone_up=arr[s1,4]*1e-6
Nsteps = arr[s1,5]*1e5
CPU = arr[s1,6]*1e6
plt.close('all')
if 0:
    fig,ax=plt.subplots(1,1)
    ax.plot(Mach,dt_sim)
    ax.plot(Mach, 1/32*1/(1+Mach))
    axbonk(ax,xlabel='Mach',ylabel='dt')
    fig.savefig('mach_dt.png')

if 1:
    fig,ax=plt.subplots(1,1)
    mus_zone_up_2=16*dt_wal/256**3
    ax.scatter(mus_zone_up, mus_zone_up_2)
    ax.plot(mus_zone_up,mus_zone_up)
    #ax.scatter(dt_wall, dt_wall/256**3)
    fig.savefig('mus_zone_up.png')
    axbonk(ax,xlabel='mus/zone/up', ylabel='mus/zone/up')
        
if 1:
    fig,ax=plt.subplots(1,1)
    #ax.scatter(Tdyn/Nsteps, dt_sim)
    #ax.scatter(dt_sim,1/(Tdyn/Nsteps/dt_sim)/10000)
    NS2 = Tdyn/dt_sim*10
    ax.scatter(Nsteps,NS2)
    ax.plot(Nsteps,Nsteps)
    axbonk(ax,xlabel='Nsteps',ylabel='Nsteps')
    fig.savefig('step_size.png')

if 1:
    fig,ax=plt.subplots(1,1)
    fac=1/56#(1e6*1/3600/64)
    print(fac)
    SU = mus_zone_up*NS2*256**3*fac
    ax.scatter(CPU,SU)
    ax.plot(CPU,CPU)
    axbonk(ax,xlabel='SU',ylabel='SU')
    fig.savefig('SU.png')
    


