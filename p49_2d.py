ef('p49_labels.py')
from mpl_toolkits.mplot3d import Axes3D
#figb, ax_RS = plt.subplots(1, 1, sharex=True)
ax_Leg=None
plt.close('all')
if 1:
    """all 4 panels"""
    figb, axes = plt.subplots(2,2,figsize=[12.8+4.8,2*4.8])
    ax_RS = axes[0][0]
    ax_MaMs = axes[0][1]
    ax_SS = axes[1][0]
    ax_Leg = axes[1][1]
    Nrows = 15
if 0:
    figb, (ax_MaMs,ax_Leg) = plt.subplots(1,2,figsize=[12.8+4.8,2*4.8])
    fig_ignore, axes = plt.subplots(2,2,figsize=[12.8+4.8,2*4.8])
    ax_RS = axes[0][0]
    #ax_MaMs = axes[0][1]
    ax_SS = axes[1][0]
    #ax_Leg = axes[1][1]
    Nrows=40
if 0:
    figb, ax_MaMs = plt.subplots(1,1,figsize=[2*4.8,2*4.8])
    fig_ignore, axes = plt.subplots(2,2,figsize=[12.8+4.8,2*4.8])
    ax_RS = axes[0][0]
    #ax_MaMs = axes[0][1]
    ax_SS = axes[1][0]
    ax_Leg = axes[1][1]
    Nrows=40

if 0:
    """just the B slope, B/E plot"""
    figb, ax_RS = plt.subplots(1,1,figsize=[2*4.8,2*4.8])
    fig_ignore, axes = plt.subplots(2,2,figsize=[12.8+4.8,2*4.8])
    #ax_RS = axes[0][0]
    ax_MaMs = axes[0][1]
    ax_SS = axes[1][0]
    ax_Leg = axes[1][1]
    Nrows=40

#ax_RS = figb.add_subplot(111, projection='3d')
#ax_RS,ax_MaMs = figb.add_subplot(1,2)

field1 = 'mach'
field2 = 'AlfMach'
labels_2d={'AlfMach':r'$M_A = \frac{V_{\rm{rms}}}{B/\sqrt{rho}}$', 'mach':r'$M_s=V_{\rm{rms}}/C_{\rm{sound}}$'}
limits = {'mach':[0.05,10.0],'AlfMach':[0.0,50]}
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
        ax_MaMs.plot( Ms, Ma, c=c,linestyle=line,label=label)
        if label is not None or beta < 0.1:
            ax_MaMs.text( 1.2*Ms[1],1.2*Ma[1], r'$\beta=%0.2f$'%beta)
        c=[0.8]*4
        ax_MaMs.plot( limits['mach'], [1.0,1.0], c=c,linestyle=":")
        ax_MaMs.plot( limits['mach'], [0.5,0.5], c=c,linestyle=":")
        ax_MaMs.plot( [1.0,1.0], limits['AlfMach'], c=c,linestyle=":")
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
title=''
if 'which_series' in dir():
    these_sims = series[which_series]
    title = titles[which_series]
    figb.suptitle(title)
#color_sim_dict = color_sim_dict_2
color_sim_dict = color_sim_dict_hsv
axis_list = 'xyz'
means = False
do_MaMs = True
modifier = 'EB_RGB_%s_'%axis_list
if means:
    modifier = 'mean_'


if 0:
    fig_names,ax_names=plt.subplots(1,1)
    names_y = 0
    names_del = 0.1
    for ncar,car_name in enumerate(these_sims):
        ax_names.text(0,ncar*names_del,car_name,color=color_sim_dict[car_name])
    outname = 'dumb_name.pdf'
    ax_names.plot([-1,-1],[0,ncar*names_del])
    ax_names.plot([1,1],[0,ncar*names_del])
    fig_names.savefig(outname)
    print(outname)

if 1:
    for ncar,car_name in enumerate(these_sims): #,'ax21']:
        print car_name
        car = taxi.taxi(car_name)
        if True and car_name in ['ac19','ac22','ac23','ac25','ac26', 'b02_512','b2_512','b20_512']:
            #car.outname = car.name + "_cmbto" #for old cmb tools
            car.outname = car.name + "_both" #new cmb tools
        car.qb_load()
        qb=car.qb

        car_list[car_name]=car
        car_outname += '%s_'%car.outname
        fs1 = field_stuff(field1, car)
        fs2 = field_stuff(field2, car)
        time = field_stuff('t_tcross', car)

        this_quan={}
        for axis in axis_list:
            Bamp, this_quana, t = scrub_eb_vs_stuff(qb,'Bamp',axis,fs1.field_name)
            Eamp, this_quana, t = scrub_eb_vs_stuff(qb,'Eamp',axis,fs1.field_name)
            Bslope, this_quana, t = scrub_eb_vs_stuff(qb,'Bslope',axis,fs1.field_name)
            Eslope, this_quana, t = scrub_eb_vs_stuff(qb,'Eslope',axis,fs1.field_name)
            Eslope2, this_quanb, t = scrub_eb_vs_stuff(qb,'Eslope',axis,fs2.field_name)

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
                #ax_MaMs.text(nominal[car_name]['mach'], nominal[car_name]['AlfMach'],'N',color=color_sim_dict.get(car.name,'k'))
                ax_MaMs.scatter(nominal[car_name]['mach'], nominal[car_name]['AlfMach'],color=color_sim_dict.get(car.name,'k')) #, label=label)
            else:
                label=None
            marker  =  {'x':'x','y':'*','z':'+'}[axis]
            #ax_SS.scatter(Bslope[sl], Eslope[sl],color=color_sim_dict.get(car.name,'k'),marker=marker) 
            ax_SS.scatter(Bslope[sl], Bslope[sl]/Eslope[sl],color=color_sim_dict.get(car.name,'k'),marker=marker) 
            the_y = Bslope[sl]
            if means:
                ax_RS.errorbar(r_bar,s_bar, xerr=r_std,yerr=s_std,color=color_sim_dict.get(car.name,'k'),marker=marker)
            else:
                r_bar = np.mean(ratios); s_bar = np.mean(the_y)
                r_std = np.std(ratios); s_std = np.std(the_y)
                ax_RS.scatter(ratios, the_y, color=color_sim_dict.get(car.name,'k'),marker=marker)
            label = None
            if axis == 'x' or len(axis_list) == 1:
                label = car.outname
            #ax_MaMs.scatter(this_quan1, this_quan2, color=color_sim_dict.get(car.name,'k'),marker=marker,label=label)
            x_bar = nar(np.mean(this_quan1)); y_bar = nar(np.mean(this_quan2))
            x_err = nar(np.std(this_quan1));  y_err = nar(np.std(this_quan2))
            if means and do_MaMs:
                ax_MaMs.errorbar(x_bar,y_bar, xerr=x_err, yerr=y_err,  color=color_sim_dict.get(car.name,'k'),marker=marker,label=label)
                pass
            elif do_MaMs:
                ax_MaMs.scatter(this_quan1, this_quan2, color=color_sim_dict.get(car.name,'k'),marker=marker,label=label)
                pass
            x1c,x2c = limits[field1]
            y1c,y2c = limits[field2]
            if max(this_quan1) > x2c:
                print "X lim too low", max(this_quan1), x2c
            if max(this_quan2) > y2c:
                print "Y lim too low", max(this_quan2), y2c
            if min(this_quan1) < x1c:
                print "X lim to high", min(this_quan1), x1c
            if min(this_quan2) < y1c:
                print "Y lim to high", min(this_quan2), y1c

            #ax_RS.plot(ratios,the_y, c=[0.5]*4)
            #ax_RS.scatter(this_quan1, this_quan2, ratios ,label=label,
            #            color=color_sim_dict[car.name])
            #for x,y in zip(this_quan1,this_quan2):
            #    ax_RS.text(x,y,axis,color=color_sim_dict[car.name])

    #ax_RS.legend(loc=0)
    x1 = 0.0; x2 = 2.0
    y1 = -10.; y2=0.
    ax_RS.set_xlabel('B/E')
    ax_RS.set_xlim(x1,x2)
    ax_RS.set_ylabel('B Slope')
    ax_RS.set_ylim(y1,y2)

    #ax_SS.plot([y1,y2],[y1,y2],c='k')
    #ax_SS.plot([y1,y2],[0.1*y1,0.1*y2],c=[0.5]*4)
    ax_SS.plot([y1,y2],[1.0,1.0], c='k')
    ax_SS.plot([y1,y2],np.ones(2)*1.5, c=[0.5]*4)
    ax_SS.plot([y1,y2],np.ones(2)*0.75, c=[0.5]*4)
    ax_SS.set_xlim(y1,y2) #; ax_SS.set_ylim(y1,y2)
    ax_SS.set_xlabel('Bslope')
    ax_SS.set_ylabel('Bslope/Eslope')

    ax_MaMs.set_xlim(x1c, x2c)
    ax_MaMs.set_ylim(y1c, y2c)
    ax_MaMs.set_xscale('log'); ax_MaMs.set_yscale('log')
    ax_MaMs.set_xlim(x1c, x2c)
    ax_MaMs.set_ylim(y1c, y2c)
    ax_MaMs.set_xlabel(labels_2d[field1])
    ax_MaMs.set_ylabel(labels_2d[field2])

    if ax_Leg is None:
        ax_MaMs.legend(loc=0)
    #ax_MaMs.legend(bbox_to_anchor=(0, 0, 1, 1), bbox_transform=gcf().transFigure)
    else:
        handles, labels = ax_MaMs.get_legend_handles_labels()
        L = []; H = []; objs=[]
        if True:
            for i in range(len(labels)/Nrows+1):
                print i
                istart = i*Nrows
                iend = (i+1)*Nrows
                dH  =  handles[istart:iend] 
                dL  =  labels[istart:iend] 
                if len(dH) :
                    L.append(dL)
                    H.append(dH)
                dO = ax_Leg.legend(dH,dL,loc=[1,2,10][i])
                objs.append(dO)
            
        for dO in  objs[:-1]:
            print 'add', dO
            plt.gca().add_artist(dO)



    ax_RS.plot([0.5]*2,[y1,y2],c=[0.5]*3)
    ax_RS.plot([0.25]*2,[y1,y2],c=[0.5]*3,linestyle=':')
    ax_RS.plot([1.0]*2,[y1,y2],c=[0.5]*3,linestyle=':')
    ax_RS.plot([x1,x2],[-2.42]*2,c=[0.5]*3)

    outname = 'p49_2d_%s_%s.png'%(modifier,car_outname)
    figb.savefig(outname)
    plt.close(figb)
    print outname

