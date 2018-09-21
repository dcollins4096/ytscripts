import matplotlib
matplotlib.use('Agg')

import yt
import matplotlib.pyplot as plt
import numpy as np
import h5py
import pdb
mython = yt.__path__[0].split('/')[-2]
print("fuu "+mython)
if 1:
#run = 'gas_plus_dm_pancake'
#run = 'unigrid'
#run = 'Adiabat_older'
    run = 'Adiabatic_Dedner'
#run = 'Adiabatic_Dedner_HackOff'
    methods=['avg','disk']
    if 'ts' not in dir() or True:
        ts = yt.load('%s/DD????/data????'%run)
        #ts = yt.load('%s/DD0006/data0006'%run)

        B_avg = np.zeros(len(ts))
        a_avg = np.zeros(len(ts))
        B_disk = np.zeros(len(ts))
        a_disk = np.zeros(len(ts))
        vel_unit = np.zeros(len(ts))
        den_unit=np.zeros(len(ts))
        mag_unit=np.zeros(len(ts))
        len_unit = np.zeros(len(ts))
        a_nath=[]
        b_nath=[]

        #for n,ds in enumerate(ts)
        n=-1
        for ds in ts:
            n=n+1
            vel_unit[n] =ds.velocity_unit.in_cgs()
            den_unit[n] =ds.mass_unit.in_cgs()/ds.length_unit.in_cgs()**3
            mag_unit[n] =ds.magnetic_unit.in_cgs()
            len_unit[n] =ds.length_unit.in_cgs()

            if 'avg' in methods:
                #B1 = np.mean(ds.all_data().quantities['WeightedAverageQuantity'](['magnetic_field_strength'],'cell_volume'))
                #B1 = np.mean(ds.all_data().quantities['WeightedAverageQuantity']([('enzo','By')],'cell_volume')).in_units('code_magnetic')
                B_avg[n] = ds.index.grids[0]['magnetic_field_strength'].mean()
                a_avg[n]= ds.scale_factor
            if 'disk' in methods:
                #B1 = np.mean(ds.all_data().quantities['WeightedAverageQuantity'](['magnetic_field_strength'],'cell_volume'))
                #B1 = np.mean(ds.all_data().quantities['WeightedAverageQuantity']([('enzo','By')],'cell_volume')).in_units('code_magnetic')
                z = ds['CosmologyCurrentRedshift']

                a = 1/(1+z) #ds.scale_factor
                a_disk[n]=a
                print(ds,z,a)
                h5_name = ds.index.grids[0].filename
                print(h5_name)
                fptr = h5py.File(h5_name,'r')
                try:
                    g = fptr['Grid%08d'%1]
                    Btot = np.sqrt(g['Bx'][:]**2+g['By'][:]**2+g['Bz'][:]**2)
                    B_disk[n] = Btot.mean()
                except:
                    raise
                finally:
                    fptr.close()
            if 'n1' in methods:
                B1 = ds.index.grids[0]['magnetic_field_strength'][0]
                print("read %s"%ds.fullpath)
                a = ds.scale_factor
                a_s[n] = a
                B_s[n] = B1
                b_nath.append(B1.mean()*a**2)
                a_nath.append(a)

    plt.clf()
#half_a = B_s*a_s**0.5
#half_a /= half_a[0]
#square_a = B_s*a_s**2
#square_a /= square_a[0]
#print("half a",half_a)
#print("a2",square_a)
    Avg_sq = B_avg*a_avg**2
    Avg_sq /= Avg_sq[0]
    disk_sq = B_disk*a_disk**2
    disk_sq /= disk_sq[0]
    plt.loglog(a_avg,Avg_sq,label='yt a^2',marker='*')
    plt.loglog(a_disk,disk_sq,label='disk a^2 ',marker='*')
    Avg_half = B_avg*a_avg**0.5
    Avg_half /= Avg_half[0]
    disk_half = B_disk*a_disk**0.5
    disk_half /= disk_half[0]
    plt.loglog(a_avg,Avg_half,label='yt a^1/2',marker='*')
    plt.loglog(a_disk,disk_half,label='disk a^1/2 ',marker='*')
    plt.legend(loc=0)
    mython = yt.__path__[0].split('/')[-2]
    plt.title('From %s'%mython)
#plt.ylim(1e-7, 1e-5)
    outdir  = '/home/dcollins4096/PigPen'
    outdir = '/Users/dcollins/RESEARCH2/EnzoProjects/F0018_cosmology_units'
    outname = '%s/%s_ad_%s_code.png'%(outdir,run, mython)
    plt.savefig(outname)
    print('out',outname)

plt.clf()
def noo(arr):
    ooo = arr/arr[0]
    print(ooo)
    return ooo
plt.plot(a_avg, noo(vel_unit),label='vel')
plt.plot(a_avg,noo(den_unit*a_avg**3), label='den a^3')
plt.plot(a_avg,noo(len_unit/a_avg), label='L a')
plt.plot(a_avg, noo( mag_unit/np.sqrt(den_unit)), label='mag/sqrt(den)')
plt.legend(loc=0)
plt.savefig('everythign_is_one.png')
