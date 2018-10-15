if 'vx' not in dir():
    ds = yt.load('/scratch/00369/tg456484/Paper42_NewAK/aq42_m9_drive0.5_512_p59_ppm/DD0040/data0040')
    ad = ds.all_data()
    vx = ad['x-velocity']
    vy, vz = ad['y-velocity'], ad['z-velocity']

vxbar = np.mean(vx, dtype='float64')
vybar = np.mean(vy, dtype='float64')
vzbar = np.mean(vz, dtype='float64')

v2 = np.mean( (vx-vxbar)**2 + (vy-vybar)**2 + (vz-vzbar)**2 )
