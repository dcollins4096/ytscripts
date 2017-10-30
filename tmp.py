if 'ef' not in dir():
    execfile('go')

reload(taxi)
junk={}
for cname in ['axb19','axb20','axb21','axb22']:
    car=taxi.taxi(cname)
    car.qb_load()
    #car.qb.plot()
    junk[car.name]=car
if 0:
    nominal = dict(zip(['ax19','ax20','ax21','ax22'],[{},{},{},{}]))
    nominal['ax19']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[1.0, 0.3,np.log10(2*(0.3/1.0)**2),11.82]))
    nominal['ax20']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[3.0, 1.0,np.log10(2*(1.0/3.0)**2),10.6347]))
    nominal['ax21']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[0.6, 0.3,np.log10(2*(0.3/0.6)**2),7.089815 ]))
    nominal['ax22']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[3.0, 0.3,np.log10(2*(0.3/3.0)**2),35.4]))

    sqrtfourpi=np.sqrt(4*np.pi)
    for cn in nominal:
        tn=nominal[cn]
        print cn, (tn['AlfMach'] - tn['mach']/(tn['field_cgs']/sqrtfourpi))/tn['AlfMach']
        print cn, "%0.2e %0.2e"%(10**(tn['logbeta']), np.pi*8/tn['field_cgs']**2)
