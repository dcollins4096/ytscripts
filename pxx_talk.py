
f1 = '/scratch1/dcollins/Paper08/B02/512/RS0060/restart0060'
f2 = '/scratch1/dcollins/Paper08/B02/512/RS0000/restart0000'
fig = plt.figure()
ax=[]
ax.append(fig.add_subplot('112'))
ax.append(fig.add_subplot('122'))
frbs = []
proj=[]
for n,f in enumerate([f1,f2]):
    ds = yt.load(f)
    proj.append(yt.ProjectionPlot(ds,0,'density'))
    ax[n]=imshow(np.log10(proj[n]._frb['density']))
    #proj.set_zlim('density',0.1,100)
    #proj.save('ppp%02d'%n)
plt.close(fig)
