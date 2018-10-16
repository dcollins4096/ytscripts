
if 0:
    dirname = '/scratch1/dcollins/Paper08/B02/512/RS%04d/restart%04d'
    ds = yt.load(dirname%(80,80))
    m,p=ds.find_max('density')
    p=p.v

if 0:
    d=[0.05*4.6]*3
    L = p-d
    R = p+d
    R[1]=1.0
    R[2] = 1.0
    box = ds.region(p,L,R)

    proj = ds.proj('density',0,data_source=box,center=p)
    print proj['density']
    pw = proj.to_pw(center=p,width = 0.1*4.6)
    print pw.save()

pw = proj.to_pw(center=p,width = 0.1*4.6)
pw.set_cmap('density','gray')
pw.annotate_magnetic_field(factor=1) #, c='b') # plot_args={'color':'b'})
print pw.save()
