#
# Spherical shell bins are not as smooth as I'd like.
# Horse around.
#


plt.close('all')
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal')
left=nar([0.,0.,0.])
right=nar([1.,1.,1.])
N = nar([32,32,32])
dx = (right-left)/N

for ix in range(N[0]+1):
    ax.plot( [ix*dx[0], ix*dx[0]], [left[1],right[1]] ,'k',linewidth=1,c=[0.5]*4)
for iy in range(N[1]+1):
    ax.plot( [left[0],right[0]] ,[iy*dx[1], iy*dx[1]], 'k',linewidth=1,c=[0.5]*4)


center = nar([0.5,0.5])
Nbins = 1024
Rinner = 0
Router = 0.25
BinSize = (Router-Rinner)/Nbins
iBin = nar(range(Nbins))
Rbin = Rinner + (0.5+iBin)*BinSize
for R in Rbin:
    circle = plt.Circle(center,R,color=[0.5,0.5,0.5],fill=False)
    ax.add_artist(circle)

outname = 'p52_stupid_grid.png'
fig.savefig(outname)
print(outname)

z,y,x=np.mgrid[left[0]+0.5*dx[0]:right[0]:dx[0], left[1]+0.5*dx[1]:right[1]:dx[1], left[2]+0.5*dx[2]:right[2]:dx[2]]
r = np.sqrt(x*x+y*y+z*z)
r_inside = r<Router
rbin = ((r[r_inside]-Rinner)/BinSize).astype('int')
rbin_flat = rbin.flatten()
r_un = np.unique(rbin)
counts = np.zeros(Nbins)
for n in r_un:
    ok = np.where( rbin_flat == n )[0].size
    counts[n]=ok
plt.clf()
#plt.plot(Rb,Nc/Nc.max()*counts.max(),label='enzo',marker='*',c='r')
plt.plot(Rbin, counts,marker='*',c='b')
outname = 'p52_counts.png'
plt.savefig(outname)


