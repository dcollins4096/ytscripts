

reload(taxi)
r441 = taxi.load('p52_r441')
r441.fields = ['density','Density_56Ni']
r441.frames=[60] #list(range(0,335,30))
r441.callbacks=['magnetic_streamlines']
r441.plot()
