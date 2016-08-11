

#
# parse random_forcing.out
#

"""
0,1     fprintf( Fptr, "%"ISYM" %9.6"FSYM" ", MetaData->CycleNumber, MetaData->Time);
2       fprintf( Fptr, " %9.3"GSYM, GlobVal[0]); // rho*v*dv  (b in the quadratic formula)
3      fprintf( Fptr, " %9.3"GSYM, GlobVal[1]); // dv*dv*rho (a in the quadratic formulat)
4      fprintf( Fptr, " %9.3"GSYM, GlobVal[2]); // v*v*rho/temp  
5      fprintf( Fptr, " %9.3"GSYM, GlobVal[3]); // v*v/temp   
6      fprintf( Fptr, " %9.3"GSYM, GlobVal[4]); // v*v*rho  (ke)
7      fprintf( Fptr, " %9.3"GSYM, GlobVal[5]); // v*v
8      fprintf( Fptr, " %9.3"GSYM, GlobVal[6]); // rho*rho
9      fprintf( Fptr, " %9.3"GSYM, GlobVal[7]); // total energy
10      fprintf( Fptr, " %9.3"GSYM, GlobVal[8]); // min density
11      fprintf( Fptr, " %9.3"GSYM, GlobVal[9]); // max density
12      fprintf( Fptr, " %9.6"FSYM, 0.50*GlobVal[4]/numberOfGridZones);   // kinetic energy
13      fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[2]/numberOfGridZones));  // mass weighted rms Mach
14      fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[3]/numberOfGridZones));  // volume weighed rms Mach
15     fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[5]/numberOfGridZones));  // rms Velocity
16     fprintf( Fptr, " %9.6"FSYM, sqrt(GlobVal[6]/numberOfGridZones));  // Density variance
17,18      fprintf( Fptr, " %9.6"FSYM" %10.5"FSYM, minDens, maxDens );                  // min/max Density
19     fprintf( Fptr, " %9.3"GSYM"\n", *norm); // rho*v*dv  (b in the quadratic formula)
"""

basedir = '/scratch/00369/tg456484/Paper42_NewAK/aq04_hydro_M9_drive_512'
basedir = '/scratch/00369/tg456484/Paper42_NewAK/aq05_hydro_m9_32'
fname = '%s/randomForcing.out'%(basedir)
sim = 'aq05'
array = np.fromfile(fname, sep=' ')
n_col = 20
lines = array.size/n_col
array.shape = (lines,n_col)
n = array[:,0]
t = array[:,1]
ke = array[:,12]
Mach_v = array[:,14]

plt.clf()
plt.plot(t,ke)
plt.plot(t,array[:,9])
plt.xlabel('time'); plt.ylabel('Kinetic Energy')
plt.savefig('%s_time_ke.pdf'%sim)
plt.clf()
plt.plot(t, Mach_v)
plt.xlabel('time'); plt.ylabel('rms mach, volume')
plt.savefig('%s_time_Mach.pdf'%sim)

plt.clf()
plt.plot(t, array[:,19])
plt.xlabel('time'); plt.ylabel('rms mach, norm')
plt.savefig('%s_time_norm.pdf'%sim)

