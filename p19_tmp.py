pos_dict={}
pos_dict['x'] = np.random.uniform(size=50)
pos_dict['y'] = np.random.uniform(size=50)
pos_dict['z'] = np.random.uniform(size=50)
denisty = 10**np.random.uniform(low=-1,high=6,size=50)
x_bar = np.mean( pos_dict['x']*density)/density.sum()
y_bar = np.mean( pos_dict['y']*density)/density.sum()
z_bar = np.mean( pos_dict['z']*density)/density.sum()
R = np.sqrt( (pos_dict['x']-x_bar)**2 + (pos_dict['y']-y_bar)**2 + (pos_dict['z']-z_bar)**2 )
r_x = pos_dict['x']-x_bar
r_y = pos_dict['y']-y_bar
r_z = pos_dict['z']-z_bar
nx = r_x/R
ny = r_y/R
nz = r_z/R
print np.sqrt(nx**2+ny**2 + nz**2)

