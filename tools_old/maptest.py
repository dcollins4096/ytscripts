#proj = yt.ProjectionPlot(ds,0,'density')
for cm in  ['BrBG']: #['Accent','BrBG','Oranges','Purples']:
  proj.set_cmap('density',cm)
  #proj.set_zlim('density',0.1,1000)
  proj.set_zlim('density',0.1,100)
  proj.save('fyap_%s'%cm)
