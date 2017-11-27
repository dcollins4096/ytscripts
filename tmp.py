if 'ef' not in dir():
    execfile('go')
if 0:
    reload(taxi)
    comp = taxi.fleet(['b02_512','b2_512','b20_512'])
    comp('car.qb_load()')
    comp('car.qb.plot(HydroMethod=6)')

if 0:

    Nrows = 3
    myarr = range(17)
    L = []; H = []
    if len(labels) > Nrows:
        for i in range(len(myarr)/Nrows+1):
            istart = i*Nrows
            iend = (i+1)*Nrows
            dummy =  myarr[istart:iend] 
            if len(dummy) :
                L.append(dummy)
            print L[-1]
    print myarr
    print L
nominal = { }
nominal['aa19']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[1  ,11.82  ,0.300 , 0.180 , -0.7]))
nominal['aa20']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[3  ,10.6347,1.000 , 0.222 , -0.7]))
nominal['aa21']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[0.6,7.089815,   0.300 , 0.500 , -0.3]))
nominal['aa22']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[3  ,35.4   ,0.300 , 0.020 , -1.7]))
nominal['ab19']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[1  ,11.82  ,0.300 , 0.180 , -0.7]))
nominal['ab22']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[3  ,35.4   ,0.300 , 0.020 , -1.7]))
nominal['ab23']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[6  ,70.8   ,0.300 , 0.005 , -2.3]))
nominal['ab24']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[1  ,330.82 ,0.011 , 0.000 , -3.6]))
nominal['ab26']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[0.1,0.07   ,5.064 , 5129.131  , 3.7]))
nominal['ac19']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[1  ,11.82  ,0.300 , 0.180 , -0.7]))
nominal['ac22']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[3  ,35.4   ,0.300 , 0.020 , -1.7]))
nominal['ac23']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[6  ,70.8   ,0.300 , 0.005 , -2.3]))
nominal['ac25']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[0.1,0.7, 0.506,  51.291, 1.7]))
nominal['ac26']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[0.1,0.07   ,5.064 , 5129.131  , 3.7]))
nominal['ax19']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[1  ,11.82  ,0.300 , 0.180 , -0.7]))
nominal['ax20']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[3  ,10.6347,1.000 , 0.222 , -0.7]))
nominal['ax21']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[0.6,7.089815,   0.300 , 0.500 , -0.3]))
nominal['ax22']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[3  ,11.82  ,0.900 , 0.180 , -0.7]))
nominal['az19']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[1  ,11.82  ,0.300 , 0.180 , -0.7]))
nominal['az20']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[3  ,10.6347,1.000 , 0.222 , -0.7]))
nominal['az21']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[0.6,7.089815,   0.300 , 0.500 , -0.3]))
nominal['az22']=  dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[3  ,35.4   ,0.300 , 0.020 , -1.7]))
nominal['b02_512']=dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[9  ,11.20998314,2.846  ,0.200 , -0.7]))
nominal['b2_512'] =dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[9  ,3.544907702,9.000  ,2.000 , 0.3]))
nominal['b20_512']=dict(zip(['mach','field_cgs','AlfMach','beta','logbeta'],[9  ,1.120998314,28.460 ,20.000, 1.3]))
axe.plot([x1,x2],[x1,x2],c='k')  
axe.plot([x1,x2],[0.1*x1,0.1*x2],c=[0.5]*4)
axe.plot([x1,x2],[1.0*x1,1.0*x2],c=[0.5]*4)
figb.savefig(outname)
