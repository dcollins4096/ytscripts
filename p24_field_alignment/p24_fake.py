axis='z'
field_horizontal = {'x':'By','y':'Bz','z':'Bx'}[axis]
field_vertical   = {'x':'Bz','y':'Bx','z':'By'}[axis]
p=1
n0=1
def _Q_local(field,data):
    """This function calculates the Stokes Parameter "Q" along an axis x, y, or z.
    Makes use of the depolarization factor "epsilon" using a power exponent.
    """
    
    epsilon = np.ones(data['density'].shape)
    n = data['density']
    B_sq = data['Bx']**2.0 + data['By']**2.0 + data['Bz']**2.0

    epsilon[ n <= n0 ] = (n)[ n <= n0 ]  
    epsilon[ n > n0 ]  = (n0**(1-p) * n**p)[ n > n0 ]   

    return ( epsilon * (((data[field_horizontal])**2.0) - ((data[field_vertical])**2.0))/B_sq )

print 'adding yt field Qz1'
yt.add_field('Qz1', units='dimensionless', function=_Q_local, force_override=True)

def _U_local(field,data):
    """Makes stokes U."""
    
    epsilon = np.ones(data['density'].shape)
    n = data['density']
    B_sq = data['Bx']**2.0 + data['By']**2.0 + data['Bz']**2.0    

    epsilon[ n <= n0 ] = (n)[ n <= n0 ]  
    epsilon[ n > n0 ]  = (n0**(1-p) * n**p)[ n > n0 ] 
    
    return  (2.0 * epsilon * ((data[field_horizontal]) * (data[field_vertical]))/B_sq)

print 'adding yt field Uz1'
yt.add_field('Uz1', units='dimensionless', function=_U_local, force_override=True)
size = [32,32,32]
Bx=np.ones(size)*0
By=np.ones(size)
By[:,:,16:] = 0.5
Bz=np.ones(size)
n =np.ones(size)
fake ={'Bx':Bx,'By':By,'Bz':Bz,'density':n}
Q = np.sum(_Q_local(None,fake), axis=2)
U = np.sum(_U_local(None,fake), axis=2)
B_hor  = np.sum(fake[field_horizontal], axis=2)
B_vert  = np.sum(fake[field_vertical], axis=2)

if 1:
    P_theta = 0.5*np.arctan2(U,Q)*180/np.pi 
    B_theta = np.arctan2(B_vert, B_hor)*180/np.pi 
#B_theta[B_theta > 90] = B_theta[B_theta > 90] - 180
#B_theta[B_theta < -90] = 180 + B_theta[B_theta < -90]

    PC_TO_PIXEL = 1
    beta = 'fake'
    res = 32
    two_pc_pix = 2.0 * PC_TO_PIXEL

    plt.clf()
    fig, ax = plt.subplots()
    plt.title("Polarization Angle Image Y %s" %beta)
    plt.imshow(P_theta, origin='lower', interpolation='nearest', cmap=plt.cm.Greys)
    plt.colorbar() 
    disk_lg = plt.Circle((res/2,res/2),two_pc_pix,color='red', fill=False)
    ax.autoscale(False)
    ax.add_patch(disk_lg)
    plt.savefig('P_theta_image_y_disk2pc_%s'%beta) 
