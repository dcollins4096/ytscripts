import numpy as np
nar=np.array
from davetools import rainbow_map
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
#nominal = {} 
#nominal['ax19']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[1.0, 0.3,np.log10(2*(0.3/1.0)**2),11.82]))
#nominal['ax20']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[3.0, 1.0,np.log10(2*(1.0/3.0)**2),10.6347]))
#nominal['ax21']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[0.6, 0.3,np.log10(2*(0.3/0.6)**2),7.089815 ]))
#nominal['ax22']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[3.0, 0.3,np.log10(2*(0.3/3.0)**2),35.4]))
#nominal['aa19']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[1.0, 0.3,np.log10(2*(0.3/1.0)**2),11.82]))
#nominal['aa20']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[3.0, 1.0,np.log10(2*(1.0/3.0)**2),10.6347]))
#nominal['aa21']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[0.6, 0.3,np.log10(2*(0.3/0.6)**2),7.089815 ]))
#nominal['aa22']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[3.0, 0.3,np.log10(2*(0.3/3.0)**2),35.4]))
#nominal['az19']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[1.0, 0.3,np.log10(2*(0.3/1.0)**2),11.82]))
#nominal['az20']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[3.0, 1.0,np.log10(2*(1.0/3.0)**2),10.6347]))
#nominal['az21']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[0.6, 0.3,np.log10(2*(0.3/0.6)**2),7.089815 ]))
#nominal['az22']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[3.0, 0.3,np.log10(2*(0.3/3.0)**2),35.4]))
#nominal['ab19']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[1.0, 0.3,np.log10(2*(0.3/1.0)**2),11.82]))
#nominal['ab22']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[3.0, 0.3,np.log10(2*(0.3/3.0)**2),35.4]))
#nominal['ab23']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[6.0, 0.3,np.log10(2*(0.3/6.0)**2),70.8]))
#nominal['abc23']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[6.0, 0.3,np.log10(2*(0.3/6.0)**2),70.8]))

#        Ma = np.sqrt(0.5*beta)*Ms
sqrtfourpi=np.sqrt(4*np.pi)
#nominal['b02_512']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[9.0, 2.8460497143030943,np.log10(0.2),3.16227786*sqrtfourpi]))
#nominal['b2_512']= dict(zip( ['mach','AlfMach','logbeta','field_cgs'],[9.0, 9.,np.log10(2),1.0*sqrtfourpi]))
#nominal['b20_512']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[9.0, 28.5,np.log10(20),0.31622854*sqrtfourpi]))

#nominal['ac19']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[1.0, 0.3,np.log10(2*(0.3/1.0)**2),11.82]))
#nominal['ac22']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[3.0, 0.3,np.log10(2*(0.3/3.0)**2),35.4]))
#nominal['ac23']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[6.0, 0.3,np.log10(2*(0.3/6.0)**2),70.8]))
#nominal['ac23']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[6.0, 0.3,np.log10(2*(0.3/6.0)**2),70.8]))
#nominal['ac25']= dict(zip(['mach','AlfMach','logbeta','field_cgs'],[0.1, 0.5,np.log10(2*(0.1/0.5)**2),0.7]))
#
#nominal['ab26'] = dict(zip(['mach','AlfMach','logbeta','field_cgs'],[0.1, 5,np.log10(2*(0.1/5)**2),0.07]))
#nominal['ac26'] = dict(zip(['mach','AlfMach','logbeta','field_cgs'],[0.1, 5,np.log10(2*(0.1/5)**2),0.07]))

Lb ='None'
Le ='None'
Mb = '*'
Me = '+'
marker_dict  =  {'x':'x','y':'*','z':'+'}
labels={}
labels['axb19']=r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('ax19', 1.0, 0.3,np.log10(2*(0.3/1.0)**2),'x')
labels['axb20']=r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('ax20', 3.0, 1.0,np.log10(2*(1.0/3.0)**2),'x')
labels['axb21']=r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('ax21', 0.6, 0.3,np.log10(2*(0.3/0.6)**2),'x')
labels['axb22']=r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('ax22', 3.0, 0.3,np.log10(2*(0.3/3.0)**2),'x')
labels['aa19'] =r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('aa19', 1.0, 0.3,np.log10(2*(0.3/1.0)**2),'x')
labels['aa20'] =r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('aa20', 3.0, 1.0,np.log10(2*(1.0/3.0)**2),'x')
labels['aa21'] =r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('aa21', 0.6, 0.3,np.log10(2*(0.3/0.6)**2),'x')
labels['aa22'] =r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('aa22', 3.0, 0.3,np.log10(2*(0.3/3.0)**2),'x')
labels['az19'] =r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('az19', 1.0, 0.3,np.log10(2*(0.3/1.0)**2),'x')
labels['az20'] =r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('az20', 3.0, 1.0,np.log10(2*(1.0/3.0)**2),'x')
labels['az21'] =r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('az21', 0.6, 0.3,np.log10(2*(0.3/0.6)**2),'x')
labels['az22'] =r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('az22', 3.0, 0.3,np.log10(2*(0.3/3.0)**2),'x')

labels['ab22'] =r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('ab22', 3.0, 0.3,np.log10(2*(0.3/3.0)**2),'x')
labels['ab19'] =r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('ab19', 1.0, 0.3,np.log10(2*(0.3/1.0)**2),'x')
labels['ab23'] =r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('ab23', 6.0, 0.3,np.log10(2*(0.3/6.0)**2),'x')
labels['b02_512'] = 'b02'
labels['b2_512' ] = 'b2'
labels['b20_512'] = 'b20'
labels['b02_512']=r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('b02', 9.0, 0.2,np.log10(2*(0.3/6.0)**2),'x')
labels['b2_512' ]=r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('b2' , 9.0, 2,  np.log10(2*(0.3/6.0)**2),'x')
labels['b20_512']=r'$\rm{%s}\ M=%0.1f\ Ma=%0.1f\ \log_{10}\beta=%0.1f\ \hat{%s}$'%('b20', 9.0, 20, np.log10(2*(0.3/6.0)**2),'x')
for lazy in ['ac19','ac22','ac23','ac25']:
    labels[lazy] = lazy
all_sims = [] # ['axb19','axb20','axb21','axb22']
all_sims += ['aa19','aa20','aa21','aa22']
all_sims += ['az19','az20','az21','az22']
all_sims += ['ab19','ab22','ab23']
all_sims += ['ac19','ac22','ac23']
all_all = np.unique(np.array(list(all_sims) + list(labels.keys())+list(nominal.keys())))
all_sims = all_all
rm = rainbow_map(len(all_sims))
color_sim_dict = dict(zip(all_sims, [rm(i) for i in range(len(all_sims))]))
color_sim_dict['aa22'] = [0.0,0.0,0.0,1.0]
for a,b in zip(['ax19','ax20','ax21','ax22'],
               ['axb19','axb20','axb21','axb22']):
    color_sim_dict[a]=color_sim_dict[b]
    labels[a]=labels[b]
color_sim_dict['b02_512']='r'
color_sim_dict['b2_512']='g'
color_sim_dict['b20_512']='b'
def scrub_eb_vs_stuff(self,eb_quantity, other_quantity, frames=[]):
    quad_frames = list(self.stuff['frames'])
    eb_frames = list(self.stuff['EBcycles'])
    if len(frames) == 0:
        frames = np.unique(nar(quad_frames + eb_frames+frames))
        frames.sort()
    ok_frames = []
    t = []
    this_quan = []
    this_eb   = []
    if other_quantity in ['brms_B']:
        other_field=(nar(self.stuff['bx2'])+nar(self.stuff['by2'])+nar(self.stuff['bz2']))/nar(self.stuff['Bx'])**2
    else:
        other_field=self.stuff[other_quantity]
    for frame in frames:
        if frame in quad_frames  and frame in eb_frames:
            ok_frames.append(frame)
            loc = np.where(self.stuff['frames']==frame)[0][0]
            loc_eb = np.where(self.stuff['EBcycles']==frame)[0][0]
            t.append(self.stuff['t'][ loc ])
            this_quan.append(other_field[loc])
            this_eb.append( self.stuff[eb_quantity][loc_eb])
    return nar(this_eb), nar(this_quan), nar(t)
class field_stuff():
    def __init__(self,field,car):
        self.car_name=car.name
        self.field=field
        self.x_scale = 'linear'
        self.logbeta=True
        self.field_name=field
        if self.field in 't_tcross':
            self.field_name = 't'
        if self.field in ['beta'] and not self.logbeta:
            self.x_scale = 'log'
    def __call__(self,array):
        out = array
        if self.field in ['beta'] and self.logbeta:
            out= np.log10(array)
        if self.field in ['t_tcross']:
            self.tcross = 0.5/nominal[self.car_name]['mach']
            out = array/self.tcross
        return out
    def x_label(self):
        out = r'$\rm{%s}$'%self.field
        if self.field in ['beta'] and self.logbeta:
            out = r'$\log_{10} \rm{%s}$'%self.field
        return out
