

#
# parse random_forcing.out
#

"""
  0,1  fprintf( Fptr, "%"ISYM" %9.6"FSYM" ", MetaData->CycleNumber, MetaData->Time);
  2  fprintf( Fptr, " %9.3"GSYM, GlobVal[0]); // rho*v*dv  (b in the quadratic formula)
  3  fprintf( Fptr, " %9.3"GSYM, GlobVal[1]); // dv*dv*rho (a in the quadratic formulat)
  4 fprintf( Fptr, " %9.3"GSYM, GlobVal[2]); // v*v*rho/temp  
  5 fprintf( Fptr, " %9.3"GSYM, GlobVal[3]); // v*v/temp   
  6 fprintf( Fptr, " %9.3"GSYM, GlobVal[4]); // v*v*rho  (ke)
  7  fprintf( Fptr, " %9.3"GSYM, GlobVal[5]); // v*v
  8        fprintf( Fptr, " %9.3"GSYM, GlobVal[6]); // rho*rho
 9  fprintf( Fptr, " %9.3"GSYM, -1e1);//GlobVal[7]); // total energy
10  fprintf( Fptr, " %9.3"GSYM, GlobVal[8]); // rho vx
11  fprintf( Fptr, " %9.3"GSYM, GlobVal[9]); // rho vy
12  fprintf( Fptr, " %9.3"GSYM, GlobVal[10]); // rho vz
13  fprintf( Fptr, " %9.3"GSYM, GlobVal[11]); // rho dvx
14  fprintf( Fptr, " %9.3"GSYM, GlobVal[12]); // rho dvy
15  fprintf( Fptr, " %9.3"GSYM, GlobVal[13]); // rho dvz
16  fprintf( Fptr, " %9.3"GSYM, GlobVal[14]); // rho
  
17  fprintf( Fptr, " %9.3"GSYM, GlobVal[GlobNum-2]); // max density
18  fprintf( Fptr, " %9.3"GSYM, GlobVal[GlobNum-1]); // max density
19  fprintf( Fptr, " %9.6"FSYM, 0.50*GlobVal[4]/numberOfGridZones);   // kinetic energy
20  fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[2]/numberOfGridZones));  // mass weighted rms Mach
21  fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[3]/numberOfGridZones));  // volume weighed rms Mach
22  fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[5]/numberOfGridZones));  // rms Velocity
23  fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[6]/numberOfGridZones));  // Density variance
245 fprintf( Fptr, " %9.6"FSYM" %10.5"FSYM, minDens, maxDens );                  // min/max Density
26  fprintf( Fptr, " %9.3"GSYM, *norm); // rho*v*dv  (b in the quadratic formula)
27  fprintf( Fptr, " %9.3"GSYM, RandomForcingVelocityOffset[0]); // rho*v*dv  (b in the quadratic formula)
28  fprintf( Fptr, " %9.3"GSYM, RandomForcingVelocityOffset[1]); // rho*v*dv  (b in the quadratic formula)
29  fprintf( Fptr, " %9.3"GSYM, RandomForcingVelocityOffset[2]); // rho*v*dv  (b in the quadratic formula)
30  fprintf( Fptr, " %9.3"GSYM, GlobVal[7]); // ke.
    fprintf( Fptr, "\n");
      fclose(Fptr);

0  1          2         3         4         5         6         7          8              9         10   11       12        13        14        15        16        17          18           19jj  
N    t        rho v dv  dvdvrho   vvrho/T   vv/t      vvrho     vv        rhorho          TE        min   max     g4 KE     mass M    vol M     Vrms      var(d)    min dens    max dens     norm
              g0        g1        g2        g3        g4        g5        g6              g7
1  0.001111   1.28e+03  1.28e+03  1.28e+03  1.28e+03  1.28e+03  1.28e+03  3.28e+04         0         1         1  0.019500  0.197484  0.197484  0.197484  1.000000  1.000000    1.00000      4.88


"""
simlist={}
simlist['aq04'] = '/scratch/00369/tg456484/Paper42_NewAK/aq04_hydro_M9_drive_512'
simlist['aq05'] = '/scratch/00369/tg456484/Paper42_NewAK/aq05_hydro_m9_32'
simlist['aq06'] = '/scratch1/dcollins/Paper42_new_turb/aq06_hydro_m9_32'
simlist['aq06b'] = '/scratch1/dcollins/Paper42_new_turb/aq06b_repeat'
simlist['aq16'] = '/scratch1/dcollins/Paper42_new_turb/aq16_ppm_m9_drive0_noamr_128'
simlist['aq16'] = '/scratch1/dcollins/Paper42_new_turb/aq15_ppm_m9_drive2_noamr_128'
simlist['aq08'] = '/scratch1/dcollins/Paper42_new_turb/aq08_hydro_m9_stationary'
simlist['aq17'] = '/scratch1/dcollins/Paper42_new_turb/aq17_ppm_m9_drive1_noamr_32'
sim='aq17'
basedir=simlist[sim]
fname = '%s/randomForcing.out'%(basedir)
array = np.fromfile(fname, sep=' ')
n_col = 32
lines = array.size/n_col
array.shape = (lines,n_col)
n = array[:,0]
t = array[:,1]
t = array[:,0]
ke = array[:,12]
Mach_v = array[:,14]
norm = array[:,29]

#plt.clf()
#plt.plot(t,array[:,19])
#plt.xlabel('time'); plt.ylabel('Kinetic Energy')
#outname = '%s_time_ke.pdf'%sim
#plt.legend(loc=0)
#plt.savefig(outname); print outname

plt.clf()
plt.plot(t,array[:,18], label='E0')
plt.xlabel('time'); plt.ylabel('max_density')
outname = '%s_dmax.pdf'%sim
#plt.ylim(-1,1)
plt.legend(loc=0)
plt.savefig(outname); print outname

plt.clf()
e0 = 0 #array[0,31]
plt.plot(t,array[:,19]-e0, label='ke')
plt.plot(t,array[:,30]-e0, label='G7')
plt.plot(t,array[:,31]-e0, label='E0')
plt.xlabel('time'); plt.ylabel('Energies')
outname = '%s_GE7.pdf'%sim
#plt.ylim(-1,1)
plt.legend(loc=0)
plt.savefig(outname); print outname

plt.clf()
tcenter = 0.5*(t[:-1]+t[1:])
dt = (t[1:]-t[:-1])
ke = array[:,19]
dedt = (ke[1:]-ke[:-1])/dt
plt.plot(tcenter,dedt, label='dedt')
target = (1.5*0.11**2-1.5*0.1**2)/0.35
plt.plot(tcenter, [target]*len(tcenter), label='target')
plt.xlabel('time (center)'); plt.ylabel('dedt')
outname = '%s_dedt.pdf'%sim
plt.legend(loc=0)
plt.savefig(outname); print outname
stat(dedt,'dedt')


plt.clf()
plt.plot(t, array[:,20],label='Mach, mass')
plt.plot(t, array[:,21],label='Mach, vol')
plt.xlabel('time'); plt.ylabel('rms mach, volume')
outname = '%s_time_Mach.pdf'%sim
plt.legend(loc=0)
plt.savefig(outname); print outname





plt.clf()
plt.xlabel('time'); #plt.ylabel('
plt.plot(t, array[:,26],label='norm')
plt.plot(t, array[:,27],label='ox')
plt.plot(t, array[:,28],label='oy')
plt.plot(t, array[:,29],label='oz')
stat( array[:,27:30], "all offset")
outname = '%s_norms.pdf'%sim
plt.legend(loc=0)
plt.savefig(outname); print outname
