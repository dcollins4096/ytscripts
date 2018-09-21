basename = "/Users/dcollins/scratch/Paper36_TracerTests/U02_2d/DD%04d/data%04d"; prefix = 'p36_U02'
#basename = "/Users/dcollins/scratch/Paper36_TracerTests/U03_Collide/DD%04d/data%04d"; prefix = 'p36_U03'

#basename = "/Users/dcollins/scratch/Paper36/t02_mhdct_square_alfv/DD%04d/data%04d"
def mnorm(field):
    return field.v - field.v.sum()/field.v.size
for n in range(11):
  ds = yt.load(basename%(n,n))
  ad = ds.all_data()
  foot = (0.0,0.0)
  oray = ds.ortho_ray(0,foot)
  #plt.plot(oray['x'], oray['density'].v-1,c='k')
  #plt.plot(oray['x'], oray['x-velocity'].v,c='r')
  plt.clf()
  by = mnorm(oray['By'])
  bz = mnorm(oray['Bz'])
  density = oray['density']
  density = density.v - density.v.sum()/density.v.size
  vx = mnorm(oray['x-velocity'])
  vy = mnorm(oray['y-velocity'])
  vz = mnorm(oray['z-velocity'])

  te = mnorm(oray['TotalEnergy'])
  #p = mnorm(oray['GasPressure'])
  plt.plot(oray['x'], density, c='k', label='d')
  plt.plot(oray['x'], vx, c='c', label = 'vx')
  #plt.plot(oray['x'], by, c='b', label='by')
  #plt.plot(oray['x'], bz, c='r', label = 'bz')
  #plt.plot(oray['x'], te, c='g', label = 'te')
  #plt.plot(oray['x'], vy, c='m', label = 'vy')
  #plt.plot(oray['x'], vz, c='y', label = 'vz')
  #plt.plot(oray['x'], p, c='c', label = 'p')
  if 1:
      px = ad['particle_position_x'].v
      pi = ad['particle_index'].v
      prange = range(px.size)
      #use this for U03
      #prange = range(0,100,10) 

      for i in prange: #range(0px.size):
          plt.text(px[i],0.01,'%s'%(int(pi[i])))
  fname = '%s_%04d'%(prefix,n)
  plt.legend(loc=2)
  plt.savefig(fname)
  print fname

