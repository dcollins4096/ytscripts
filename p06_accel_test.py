

Edges=[-4.0,4.0]
Nz =16
dx = (Edges[1]-Edges[0])/Nz
z=np.arange(Edges[0]-5*dx,Edges[1]+5*dx,dx, dtype='float')
z_half=np.arange(Edges[0]-5*dx,Edges[1]+5*dx,dx, dtype='float')+0.5*dx
center = 0.0

