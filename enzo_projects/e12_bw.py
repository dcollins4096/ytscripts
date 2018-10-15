
basedir='/Users/dcollins/scratch/EnzoProjects/E12/'
sim1 = basedir + 'd03b_ct_bw_noCT/DD%04d/data%04d'%(11,11)
sim2 = basedir + 'd03_ct_bw/DD%04d/data%04d'%(11,11)
fields=['density','pressure','TotalEnergy','x-velocity', 'By' ]#,'x-velocity']
ds1=yt.load(sim1)
ds2=yt.load(sim2)
y=tube.tube([ds1,ds2],  delta=False, fields=fields, renorm=False,filename = 'e12_tube.pdf',
            legend=True, labels=['noct','ct'])
