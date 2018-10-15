ds1 = yt.load('/scratch1/dcollins/Paper48_DrivingFiller/u05/DD0000_fresh/data0000')
ds2 = yt.load('/scratch1/dcollins/Paper48_DrivingFiller/Bred8mhd/DD0001/data0001')

for n,ds in enumerate([ds1,ds2]):
    plot = yt.ProjectionPlot(ds,0,'density')
    plot.save("p48_test_%d"%n)
