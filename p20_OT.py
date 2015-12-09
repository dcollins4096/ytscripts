if 'ef' not in dir():
    execfile('go')
if 'frame' not in dir():
    frame = 0
basedir = '/scratch1/dcollins/Paper20/a05_OT'
basedir = '/scratch1/dcollins/Paper20/a08_angle_test'
ef('p20_gradN.py')
prefix = 'a08_angle'
line_of_sight = 2
setname = '%s/DD%04d/data%04d'%(basedir,frame,frame)
ds = yt.load(setname)
proj_den = ds.proj('density',line_of_sight)
size = [128]*2
frb_den = proj_den.to_frb(1,size)
proj_mag= ds.proj('Bx',line_of_sight, weight_field = 'cell_mass') 
frb_mag = proj_mag.to_frb(1,size)
#make_plotses(frb_den,frb_mag,'test_OT')
plt.clf()
plt.imshow(frb_den['density'].v,origin='lower',interpolation='nearest',cmap='gray')
plt.colorbar()
nx,ny = frb_den['density'].shape
Y,X = np.mgrid[0:nx:1, 0:ny:1]
B1='Bx'
B2='By'
plt.streamplot(X,Y, frb_mag[B1].v, frb_mag[B2].v,color='b')
plt.xlim([0,nx])
plt.ylim([0,ny])
plt.savefig('%s_dq.png'%prefix)
Q = frb_den['Q']
U = frb_den['U']
formt = '%0.2f'
stat(frb_den['By'],"By", format=formt)
stat(frb_den['Bx'],"Bx", format=formt)
stat(Q, "Q!", format=formt)
stat(U, "U!", format=formt)

theta_derp = 0.5*np.arctan(U/Q)*180/np.pi
stat(theta_derp, "Theta Stokes plain arctan", format=formt)
theta = 0.5*np.arctan2(U,Q)*180/np.pi
stat(theta, "Theta Stokes 2",format=formt)
theta_b=np.arctan(frb_den['By']/frb_den['Bx'])*180/np.pi
stat(theta_b, "Theta B",format=formt)

plt.clf()
plt.imshow(theta, origin='lower', interpolation='nearest',cmap='hsv')
plt.colorbar()
plt.savefig('%s_stokes.png'%prefix)
plt.clf()
plt.imshow(theta_b, origin='lower', interpolation='nearest',cmap='hsv')
plt.colorbar()
plt.savefig('%s_theta_direct.png'%prefix)
