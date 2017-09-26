

Edges=[-4.0,4.0]
Nz =16
dx = (Edges[1]-Edges[0])/Nz
z=np.arange(Edges[0]-5*dx,Edges[1]+5*dx,dx, dtype='float')
z_half=np.arange(Edges[0]-5*dx,Edges[1]+5*dx,dx, dtype='float')+0.5*dx
z_mhalf=np.arange(Edges[0]-5*dx,Edges[1]+5*dx,dx, dtype='float')-0.5*dx
GravConst =6.67428e-8
az = -1200*GravConst*np.sign(z_mhalf)
Dims2 = z.size
k=np.arange(Dims2)
center = 0.0
BoxCenter = 0.0
k_center = (BoxCenter-z.min())/dx
top_bool = (k + k_center) < Dims2
shift_top = -1
bottom_bool = (k + k_center) >= Dims2
top_k = ((k + k_center)[top]).astype('int')
bottom_k = ((Dims2 - k - 1)[bottom]).astype('int')

shift_bottom = 1
b_args = np.argsort(b[:,0])
b_srt = b[b_args]
az_code = b_srt[:,1]
integrated_z = np.zeros_like(az)
integrated_z[top_k] = np.cumsum(-1*shift_top*az[top_k]*(dx))
integrated_z[bottom_k] = np.cumsum(-1*shift_bottom*az[bottom_k]*(dx))+integrated_z[13]
if 0:
    #print np.abs(az_code - az).sum()/np.mean(np.abs(az_code))
    #print integrated_z[top_k]- b_srt[:,2][top_k]
    print integrated_z[top_k]
    print b_srt[:,2][top_k]
    print np.abs(integrated_z[top_k]- b_srt[:,2][top_k]) < 1e-16
    print integrated_z[bottom_k]
    print b_srt[:,2][bottom_k]
    print np.abs(integrated_z[bottom_k]- b_srt[:,2][bottom_k]) < 1e-16
    #print b_srt[:,6]

max_int_code = b_srt[:,3]
max_integrated_accel= np.abs(integrated_z).max()
Gamma = 1.666

StratifiedBoxUniformGasEnergy = 0.016204238299692
StratifiedBoxUniformDensity = 0.59786381703727 
P0 = StratifiedBoxUniformGasEnergy*StratifiedBoxUniformDensity*(Gamma-1.0)
K0 = P0/(StratifiedBoxUniformDensity**(Gamma))
A0 = ( (Gamma-1.0)/(K0*Gamma) )**(Gamma-1.0)
CentralDensity = A0*(max_integrated_accel**(1./(Gamma-1.0))) + StratifiedBoxUniformGasEnergy;
P01= 0.0064521598895156 
A01= 12.42926005292 
K01= 0.0090884120498036

#print "rhoc relerr",1-CentralDensity / b_srt[:,5][0]
density = np.zeros_like(az)
density[top_k] = np.sign(integrated_z[top_k]) * A0*np.abs(integrated_z[top_k])**(1./(Gamma-1)) + CentralDensity
density[bottom_k] = np.sign(integrated_z[bottom_k]) * A0*np.abs(integrated_z[bottom_k])**(1./(Gamma-1)) + CentralDensity
p = K01*density**Gamma
#print np.abs(1-density / b_srt[:,7]) <1e-12

fz = density * az
gradp = -(p[1:]-p[:-1])/dx

