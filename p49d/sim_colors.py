#color is Alfven
#line is Sonic
#marker is Sonic
from GL import *

#these need to stay in this order.
simlist=nar(["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"])
color_list=nar(['r','g','b','r','g','b','r','g','b','r','g','b'])
line_list=nar(['-','-','-','-.','-.','-.','--','--','--',':',':',':'])
marker_list = nar(['.','.','.','s','s','s','^','^','^','*','*','*'])


sim_ms = nar(['half','1','2','3'])
sim_ma = nar(['half','1','2'])
plot_order=[]
for ma in sim_ma:
    for ms in sim_ms:
        sim="%s_%s"%(ms,ma)
        plot_order.append(sim)

color=dict(zip(simlist,color_list))
linestyle = dict(zip(simlist,line_list))
marker = dict(zip(simlist,marker_list))

def vals_from_sim(sim):
    ms,ma = sim.split("_")
    if ms == 'half':
        ms = 0.5
    if ma == 'half':
        ma = 0.5
    ms=float(ms)
    ma=float(ma)
    return ms,ma

ms_list=[]
ma_list=[]
for sim in simlist:
    ms,ma = vals_from_sim(sim)
    ms_list.append( ms)
    ma_list.append(ma)
ms_list=nar(ms_list)
ma_list=nar(ma_list)

Means = {'1_1':  nar([3.55, -2.52619106e-17,9.53081472e-18]),
      '1_half':  nar([7.09, -7.90451146e-18,-3.01035509e-18]),
         '2_2':  nar([3.54, 8.08340482e-17, 5.56873341e-17]),
         '3_1':  nar([10.6, 3.73914224e-17, -2.14807555e-17]),
      '3_half':  nar([21.3, 1.51618898e-17, 7.8875708e-18]),
      'half_2':  nar([1.77, 8.93619759e-19, 3.89635156e-19]),
         '1_2':  nar([1.78, 6.06340065e-17, -1.99222149e-18]),
         '2_1':  nar([7.09, 2.53215417e-16, -3.00188477e-17]),
      '2_half':  nar([14.2, -6.57297567e-19,1.95834017e-18]),
         '3_2':  nar([5.32, -1.48427277e-16,8.21283146e-17]),
      'half_1':  nar([1.77, 8.93619759e-19, 3.89635156e-19]),
   'half_half':  nar([3.54, 1.3471212e-17,  2.01932655e-18])}

