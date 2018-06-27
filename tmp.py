my_ts = tsA01 #tsy701
my_stuff = stuff
directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/rA01_rb96_110_f-/%s'
#directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/y701_rb96_fft_f-_play/%s'
off_disk = {}
names={}
n2={}
n2['hx']='Bx'
n2['hy']='By'
n2['hz']='Bz'
n2['d']='density'
n2['vx']='x-velocity'
n2['vy']='y-velocity'
n2['vz']='z-velocity'
names['hx']='Bx_16.h5'
names['hy']='By_16.h5'
names['hz']='Bz_16.h5'
names['d']='density_16.h5'
names['vx']='x-velocity_16.h5'
names['vy']='y-velocity_16.h5'
names['vz']='z-velocity_16.h5'
#['d','density'],['vx','x-velocity'],['hx','Bx'],['hy','By'],['hz','Bz'],
#                ['vz','z-velocity'],['vy','y-velocity'], ['p','GasPressure']
for field in ['hx']: #['d','vx','vy','vz','hx','hy','hz']:# field_list:
    #print( field, np.abs(my_ts.temp_means[field])-my_stuff['means'][field])

    lc = my_stuff['cubes'][field]
    oc = my_ts.temp_cubes[field]
    occ = my_ts.cubes[n2[field]]

    df = h5py.File(directory%names[field])
    dc_all = df[names[field]][:]
    dcc = dc_all.swapaxes(0,2)
    dc = dc_all.swapaxes(0,2)[:16,:16,:16]
    df.close()

    ff = h5py.File(directory%('DD0000/data0000.cpu0000'))
    fc_all = ff['Grid00000001']['BxF'][:]
    fcc = fc_all.swapaxes(0,2)
    fc = fc_all.swapaxes(0,2)[:16,:16,:16]
    ff.close()
    gf = h5py.File(directory%('DD0000/data0000.cpu0000'))
    gc_all = gf['Grid00000001']['Bx'][:]
    gcc = gc_all.swapaxes(0,2)
    gc = gc_all.swapaxes(0,2)[:16,:16,:16]
    ff.close()
    
    bar = 0.# np.mean(oc)
    print( field, dif)
    print("lc  ",lc[:5,8,8]-bar)
    print("occ ",occ[:5,8,8]-bar)
    print("dcc ",dcc[:5,8,8]-bar)
    print( "dc - oc %0.2e"%(np.abs(dc-oc).sum()))
    print( "dc - lc %0.2e"%(np.abs(dc-lc).sum()))
    print( "dcc - occ %0.2e"%(np.abs(dcc-occ).sum()))
