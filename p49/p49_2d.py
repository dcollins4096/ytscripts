from go import *
#from p49_stuff import *
from p49_labels import *
from mpl_toolkits.mplot3d import Axes3D
these_sims = ['ax19','ax22','az19','az22']
#figb, axb = plt.subplots(1, 1, sharex=True)
axd=None
plt.close('all')
if 0:
    figb, axes = plt.subplots(2,2,figsize=[12.8+4.8,2*4.8])
    axb = axes[0][0]
    axc = axes[0][1]
    axe = axes[1][0]
    axd = axes[1][1]
    figc = figd = fige = None
else:
    figb, axb = plt.subplots(1)
    figc, axc = plt.subplots(1)
    figd, axd = plt.subplots(1)
    fige, axe = plt.subplots(1)
counter=0

#axb = figb.add_subplot(111, projection='3d')
#axb,axc = figb.add_subplot(1,2)

field1 = 'mach'
field2 = 'AlfMach'
limits = {'mach':[0.05,10.0],'AlfMach':[0.0,50]}
car_list={}
car_outname = ""
these_sims = ['aa19','aa20','aa21','aa22']
these_sims += ['az19','az20','az21','az22']
these_sims += ['ax19','ax20','ax21','ax22']
#these_sims = ['ax19','ax20','ax21','ax22']
#these_sims = ['aa19','aa20','aa21','aa22']
#these_sims = ['az19','az20','az21','az22']
#these_sims = ['ax22','aa22','az22']
#these_sims = ['ax19','aa19','az19']
#these_sims = ['ax21','aa21','az21']
#these_sims = ['ax20','aa20','az20']
these_sims = ['ab19','ab22','ab23']
these_sims += ['aa19','aa20','aa21','aa22']
#these_sims = ['ax21']# 'ax19','ax20',,'ax22']
#these_sims =['ab19']
#these_sims = ['ab19','aa19','ax19','az19']
these_sims = ['ab19','aa19','ax19','az19']
these_sims = ['ab22','aa22','ax22','az22']
these_sims = ['b02_512','b2_512','b20_512']
these_sims =  ['ac19','ac22'] #'ac23','ac25']
these_sims += ['ax19','ax22'] #'ac23','ac25']
#these_sims = all_sims
these_sims = all_sims
#these_sims = ['ab26','ac26','ac23','ac25', 'ac22']

#b  Supersonic, self-gravitating b02_512 b2_512 b20_512
#ax Supersonic, trans&sub-alfvenic, Isothermal, crashed axb19 axb20 axb21 axb22 
#  Supersonic, trans&sub, 512, gamm=5/3 aa19 aa20 aa21 aa22
#  Supersonic, trans&sub, 256, gamm=5/3 ab19 ab22 ab23 ab24 ab26
#  Supersonic, trans&sub, 512, gamma=1.001, ac19 ac22 ac23 ac25 ac26
#  ayb20 ayb22
#  az19 az20 az21 az22


these_sims = ['b02_512','b2_512','b20_512']
these_sims = ['ax22', 'aa22','ab22','ac22']
these_sims = all_sims
these_sims=['aa19', 'aa20', 'aa21', 'aa22', 'ab19', 'ab22', 'ab23', 'ab24',
       'ab26', 'ac19', 'ac22', 'ac23', 'ac25', 'ac26', 'ax19', 'ax20',
              'ax21', 'ax22', 'az19', 'az20',
                     'az21', 'az22', 'b02_512', 'b20_512', 'b2_512']
#these_sims = ['aa19','aa20','aa21','aa22']
#these_sims = ['ac23']
TheGoodSims=['ab25','ab27','ab28','ab29', 'ac19', 'ac22', 'ac23', 'ac25', 'ac26', 'b02_512', 'b20_512', 'b2_512']
TheGoodSims=['ab25','ab28','ab29', 'ac19', 'ac22', 'ac23', 'ac25', 'ac26', 'b02_512', 'b20_512', 'b2_512']
these_sims=TheGoodSims
if 1:
    color_sim_dict={}
    rm = rainbow_map(len(these_sims))
    for n,car_name in enumerate(these_sims): #,'ax21']:
        color_sim_dict[car_name]=rm(n)


#color_sim_dict = dict(zip(['aa19','aa20','aa21','aa22'],['r','g','b','k']))
axis_list = 'xyz'
means = True
modifier = 'all_'
if means:
    modifier = 'mean_'
modifier = 'take11_'+axis_list
clear_plots=False
if 1:
    for ncar,car_name in enumerate(these_sims): #,'ax21']:
        print(car_name)
        if clear_plots:
            axb.clear()
            axc.clear()
            axe.clear()
            axd.clear()
        car = taxi.taxi(car_name)
        car.qb_load('p49_data_h5/quan_box_%s.h5'%car.outname)
        qb=car.qb

        car_list[car_name]=car
        car_outname += '%s_'%car.outname
        fs1 = field_stuff(field1, car)
        fs2 = field_stuff(field2, car)
        time = field_stuff('t_tcross', car)

        this_quan={}
        for axis in axis_list:
            Bamp, this_quana, t = scrub_eb_vs_stuff(qb,'Bamp_%s'%axis,fs1.field_name)
            Eamp, this_quana, t = scrub_eb_vs_stuff(qb,'Eamp_%s'%axis,fs1.field_name)
            Bslope, this_quana, t = scrub_eb_vs_stuff(qb,'Bslope_%s'%axis,fs1.field_name)
            Eslope, this_quana, t = scrub_eb_vs_stuff(qb,'Eslope_%s'%axis,fs1.field_name)
            Eslope2, this_quanb, t = scrub_eb_vs_stuff(qb,'Eslope_%s'%axis,fs2.field_name)

            tcross = time(t)

            this_quan1 = fs1(this_quana)
            this_quan2 = fs2(this_quanb)
            sl=np.arange(len(Bamp))
            mask1 = np.logical_and(np.isnan(Bamp) == False,np.isnan(this_quan1) == False)
            mask1 = np.logical_and( mask1, np.isnan(this_quan2) == False)
            if max(tcross) > 1:
                mask1 = np.logical_and(mask1, tcross>0.9)
            sl = sl[mask1]
            this_quan1 = this_quan1[sl]
            this_quan2 = this_quan2[sl]
            #sl=slice(None)
            ratios = Bamp[sl]/Eamp[sl]
            ratios[ ratios > 1.5 ] = 1.5
            if axis == 'x' or len(axis_list) == 1:
                label = car.name
                #axc.text(nominal[car_name]['mach'], nominal[car_name]['AlfMach'],'N',color=color_sim_dict.get(car.name,'k'))
            else:
                label=None
            marker  =  {'x':'x','y':'*','z':'+'}[axis]
            #axe.scatter(Bslope[sl], Eslope[sl],color=color_sim_dict.get(car.name,'k'),marker=marker) 
            axe.scatter(Bslope[sl], Bslope[sl]/Eslope[sl],color=color_sim_dict.get(car.name,'k'),marker=marker) 
            the_y = Bslope[sl]
            if means:
                axb.scatter(ratios, the_y, color=color_sim_dict.get(car.name,'k'),marker=marker)
                pass
            else:
                r_bar = np.mean(ratios); s_bar = np.mean(the_y)
                r_std = np.std(ratios); s_std = np.std(the_y)
                axb.scatter(ratios, the_y, color=color_sim_dict.get(car.name,'k'),marker=marker)
                axb.errorbar(r_bar,s_bar, xerr=r_std,yerr=s_std,color=color_sim_dict.get(car.name,'k'),marker=marker)
            label = None
            if axis == 'x' or len(axis_list) == 1:
                label = car.outname
            #axc.scatter(this_quan1, this_quan2, color=color_sim_dict.get(car.name,'k'),marker=marker,label=label)
            x_bar = nar(np.mean(this_quan1)); y_bar = nar(np.mean(this_quan2))
            x_err = nar(np.std(this_quan1));  y_err = nar(np.std(this_quan2))
            if means:
                axc.scatter(this_quan1, this_quan2, color=color_sim_dict.get(car.name,'k'),marker=marker,label=label)
            else:
                #axc.scatter(x_bar, y_bar, color=color_sim_dict.get(car.name,'k'),marker=marker,label=label)
                axc.errorbar(x_bar,y_bar, xerr=x_err, yerr=y_err,  color=color_sim_dict.get(car.name,'k'),marker=marker,label=label)
            x1c,x2c = limits[field1]
            y1c,y2c = limits[field2]
            if max(this_quan1) > x2c:
                print("X lim too low", max(this_quan1), x2c)
            if max(this_quan2) > y2c:
                print("Y lim too low", max(this_quan2), y2c)
            if min(this_quan1) < x1c:
                print("X lim to high", min(this_quan1), x1c)
            if min(this_quan2) < y1c:
                print("Y lim to high", min(this_quan2), y1c)

            #axb.plot(ratios,the_y, c=[0.5]*4)
            #axb.scatter(this_quan1, this_quan2, ratios ,label=label,
            #            color=color_sim_dict[car.name])
            #for x,y in zip(this_quan1,this_quan2):
            #    axb.text(x,y,axis,color=color_sim_dict[car.name])

        #axb.legend(loc=0)
        ###
        if field1 == 'mach' and field2 == 'AlfMach':
            for nb, beta in enumerate([0.01,0.1, 0.5, 1.0,  2,  10, 100, 1000]):
                Ms = nar(limits['mach'])

# beta = P/B^2/2
#      = 2 cs^2 rho/B^2
#      = 2 Ma^2 / Ms^2
# Ms sqrt(beta/2) = Ma
                Ma = np.sqrt(0.5*beta)*Ms
                line='-'
                c=[0.5]*4
                label = None
                beta_string = "%0.2f"%beta
                if np.abs(beta - 1.0) < 1e-5:
                    c='k'
                    label = r'$\beta=%s$'%beta_string
                if np.abs(beta - 0.1) < 1e-5: 
                    line = '--'
                    label = r'$\beta=%s$'%beta_string
                if np.abs(beta - 10) < 1e-5:
                    line = ":"
                    label = r'$\beta=%s$'%beta_string
                axc.plot( Ms, Ma, c=c,linestyle=line,label=label)
        ###
        x1 = 0.0; x2 = 2.0
        y1 = -10.; y2=0.
        axb.set_xlabel('B/E')
        axb.set_xlim(x1,x2)
        axb.set_ylabel('B Slope')
        axb.set_ylim(y1,y2)

        #axe.plot([y1,y2],[y1,y2],c='k')
        #axe.plot([y1,y2],[0.1*y1,0.1*y2],c=[0.5]*4)
        axe.plot([y1,y2],[1.0,1.0], c='k')
        axe.plot([y1,y2],np.ones(2)*1.5, c=[0.5]*4)
        axe.plot([y1,y2],np.ones(2)*0.75, c=[0.5]*4)
        axe.set_xlim(y1,y2) #; axe.set_ylim(y1,y2)
        axe.set_xlabel('Bslope')
        axe.set_ylabel('Bslope/Eslope')

        axc.set_xlim(x1c, x2c)
        axc.set_ylim(y1c, y2c)
        axc.set_xscale('log'); axc.set_yscale('log')
        axc.set_xlim(x1c, x2c)
        axc.set_ylim(y1c, y2c)
        axc.set_xlabel(field1)
        axc.set_ylabel(field2)

        if axd is None:
            axc.legend(loc=0)
        #axc.legend(bbox_to_anchor=(0, 0, 1, 1), bbox_transform=gcf().transFigure)
        elif None:
            handles, labels = axc.get_legend_handles_labels()
            Nrows = 15
            L = []; H = []; objs=[]
            if True:
                for i in range(len(labels)//Nrows+1):
                    print(i)
                    istart = i*Nrows
                    iend = (i+1)*Nrows
                    dH  =  handles[istart:iend] 
                    dL  =  labels[istart:iend] 
                    if len(dH) :
                        L.append(dL)
                        H.append(dH)
                    dO = axd.legend(dH,dL,loc=[1,2,10][i])
                    objs.append(dO)
                
            for dO in objs[:-1]:
                print('add', dO)
                plt.gca().add_artist(dO)



        axb.plot([0.5]*2,[y1,y2],c=[0.5]*3)
        axb.plot([0.25]*2,[y1,y2],c=[0.5]*3,linestyle=':')
        axb.plot([1.0]*2,[y1,y2],c=[0.5]*3,linestyle=':')
        axb.plot([x1,x2],[-2.42]*2,c=[0.5]*3)

        counter += 1
        if clear_plots:
            if figc is None:
                pdb.set_trace()
                outname = 'p49_plots_temp/p49_2d_%s_%s.png'%(modifier,car_name)
                #outname = 'p49_plots_temp/p49_stuff_%d.png'%counter
                figb.savefig(outname)
                print(outname)
            else:
                outname_base = 'p49_plots_temp/p49_ind'
                figb.savefig(outname_base+car_name +"_b" +".png")
                figc.savefig(outname_base+car_name +"_c" +".png")
                #figd.savefig(outname_base+"_d" +car_name +".png")
                #fige.savefig(outname_base+"_e" +car_name +".png")
                print(outname_base)
if not clear_plots:
    if figc is None:
        figb.savefig(outname)
        outname = 'p49_plots_temp/p49_stuff_%d.png'%counter
        print(outname)
    else:
        outname_base = 'p49_plots_temp/p49_ind_%d'%counter
        figb.savefig(outname_base + "_b.png")
        figc.savefig(outname_base + "_c.png")
        #figd.savefig(outname_base + "_d.png")
        #fige.savefig(outname_base + "_d.png")
        print(outname_base)
#plt.close(figb)

