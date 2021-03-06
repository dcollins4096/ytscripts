if 'ef' not in dir():
    execfile('go')
ef('turb_quan.py')
ef('p49_stuff.py')
qb_list = []
car_list=[]
car_outname = ''
plot_format='png'
out_format= 'p49_%s_%s_corr_tcut_%s.%s'
if 1: #'horz' not in dir():
    horz=[]
    rats=[]
    Baaa=[]
    Eaaa=[]
    Esss=[]
    Bsss=[]
    axxx=[]
    marker = []
    c=[]
    plots_to_annotate=[]
    annotations=[]

these_sims = ['aa19','aa20','aa21','aa22']
these_sims += ['az19','az20','az21','az22']
these_sims += ['ax19','ax20','ax21','ax22']
#these_sims = ['ax22','aa22','az22']
#these_sims = ['ax19','aa19','az19']
#these_sims = ['ax21','aa21','az21']
these_sims = ['ax22','aa22','az22']

#these_sims = ['ab19','aa19','ax19','az19']
#these_sims = ['b02','b2','b20']
these_sims = all_sims

#field= 't'
#field= 't_tcross'
#field = 'AlfMach'
#field = 'beta'
#field = 'mach'
field = 'brms_B'
axis_list = 'xyz'
modifier = 'take5_'+axis_list
these_sims = ['b02_512','b2_512','b20_512']
these_sims=['aa19', 'aa20', 'aa21', 'aa22', 'ab19', 'ab22', 'ab23', 'ab24',
       'ab26', 'ac19', 'ac22', 'ac23', 'ac25', 'ac26', 'ax19', 'ax20',
              'ax21', 'ax22', 'az19', 'az20',
                     'az21', 'az22', 'b02_512', 'b20_512', 'b2_512']

if 0:
    fig, ax90210 = plt.subplots(3, 1, sharex=True)
    (Rx_amp, Rx_slope, Rx_ratio)  = ax90210
    fig = [fig]
    fig_names = [modifier+'all3']
else:
    fig_amp, Rx_amp = plt.subplots(1, 1)
    fig_slope, Rx_slope = plt.subplots(1, 1)
    fig_ratio, Rx_ratio = plt.subplots(1, 1)
    figs = [fig_amp,fig_slope,fig_ratio]
    fig_names = [modifier+'amp',modifier+'slope',modifier+'rato']
    ax90210 = [Rx_amp,Rx_slope,Rx_ratio]



if 1:
    for ncar,car_name in enumerate(these_sims): #,'ax21']:
        car = taxi.taxi(car_name)
        qb = quan_box(car, plot_format=plot_format)
        qb.load()
        car_list.append(car)
        qb_list.append(qb)
        car_outname += '%s_'%car.outname
        fs = field_stuff(field, car)
        time = field_stuff('t_tcross', car)

        for axis in axis_list:
            Bamp, this_quana, t = scrub_eb_vs_stuff(qb,'Bamp',axis,fs.field_name)
            Eamp, this_quana, t = scrub_eb_vs_stuff(qb,'Eamp',axis,fs.field_name)
            Bslope, this_quana, t = scrub_eb_vs_stuff(qb,'Bslope',axis,fs.field_name)
            Eslope, this_quana, t = scrub_eb_vs_stuff(qb,'Eslope',axis,fs.field_name)
            this_quan = fs(this_quana)
            tcross = time(t)
            if field in ['t', 't_tcross']:
                to_plot=fs(t)
                scatter=False
            else:
                to_plot=fs(this_quana)
                scatter=True

            sl=np.arange(len(Bamp))
            mask1 = np.logical_and(np.isnan(Bamp) == False,np.isnan(this_quan) == False)
            if max(tcross) > 1:
                mask1 = np.logical_and(mask1, tcross>0.9)
            sl = sl[mask1]
            tcross=tcross[sl]

            #sl=slice(None)
            ratios = Bamp[sl]/Eamp[sl]
            ratios[ ratios > 1.5] = 1.5

            Npoints = len(to_plot[sl])
            horz += to_plot[sl].tolist()
            rats += (ratios).tolist()
            Baaa += Bamp[sl].tolist()
            Eaaa += Eamp[sl].tolist()
            Bsss += Bslope[sl].tolist()
            Esss += Eslope[sl].tolist()
            axxx += [axis] * Npoints
            c     = color_sim_dict[car_name]# ]*Npoints
            #c='k'
            Lb ='None'
            Le ='None'
            Mb = '*'
            Me = '+'
            marker  =  {'x':'x','y':'*','z':'+'}[axis]
            this_label=None
            if axis == 'x' or len(axis_list) == 1:
                this_label = labels[car_name]

            Rx_amp.scatter(to_plot[sl],Bamp[sl], c=c, marker=marker, label=this_label)
            if len(figs) > 1:
                more_label = this_label
            else:
                more_label=None
            Rx_amp.scatter(to_plot[sl],Eamp[sl], c=c, marker=marker, label=more_label)
            Rx_slope.scatter(to_plot[sl],Bslope[sl], c=c, marker=marker, label=more_label)
            Rx_slope.scatter(to_plot[sl],Eslope[sl], c=c, marker=marker)

            Rx_ratio.scatter(to_plot[sl],ratios, c=c, marker=marker, label=more_label)
            if 1:
                Rx_amp.plot(to_plot[sl],Bamp[sl],     c=[0.5]*4)# marker=marker)
                Rx_amp.plot(to_plot[sl],Eamp[sl],     c=[0.5]*4)# marker=marker)
                Rx_slope.plot(to_plot[sl],Bslope[sl], c=[0.5]*4)# marker=marker)
                Rx_slope.plot(to_plot[sl],Eslope[sl], c=[0.5]*4)# marker=marker)
                Rx_ratio.plot(to_plot[sl],ratios,     c=[0.5]*4)# marker=marker)

Baaa=nar(Baaa)
Eaaa=nar(Eaaa)
Bsss=nar(Bsss)
Esss=nar(Esss)
rats=nar(rats)
horz=nar(horz)


scatter=True
ext_amp     = collect_extrema(Baaa, None)
ext_amp     = collect_extrema(Eaaa, ext_amp)
ext_slope   = collect_extrema(Bsss, None)
ext_slope   = collect_extrema(Esss, ext_slope)
ext_ratio   = collect_extrema(rats, None)
ext_horz    = collect_extrema(horz, None)

if 1:
    #nominal values
    Rx_slope.plot( [horz.min(),horz.max()],[-2.42]*2,c=[0.5]*3, label='target')
    Rx_ratio.plot( [horz.min(),horz.max()],[0.5]*2,c=[0.5]*3) #B = 0.5 E.  More e.
    Rx_ratio.plot( [horz.min(),horz.max()],[1.0]*2,c=[0.5]*3, linestyle=':') #B = 0.5 E.  More e.
    Rx_ratio.plot( [horz.min(),horz.max()],[0.25]*2,c=[0.5]*3, linestyle=':') #B = 0.5 E.  More e.
for ax in ax90210:
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
Rx_amp.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#if field in nominal['ax19'].keys()
#    for sim in these_sims:


def bump2(extents, log=False):
    out = nar([0.0,0.0])
    if log:
        out[0] =10**(np.floor( np.log10(extents[0])))
        out[1] =10**(np.ceil(  np.log10(extents[1])))
    else:
        out[0] = np.floor(extents[0])
        out[1] = np.ceil(extents[1])
        #print out, extents, factor
    return out
def bump1(extents, factor, log=False):
    out = nar([0.0,0.0])
    if log:
        out[0] =np.log(extents[0])- factor*np.abs(np.log(extents[0]))
        out[1] =np.log(extents[1])+ factor*np.abs(np.log(extents[1]))
        out[0] = np.exp(out[0])
        out[1] = np.exp(out[1])
    else:
        out[0] = extents[0] - factor*np.abs(extents[0])
        out[1] = extents[1] + factor*np.abs(extents[1])
        #print out, extents, factor
    return out
bump = bump2


if 0:
    #Rx_amp.set_ylim(  0, 1.0)
    Rx_slope.set_ylim( min([-8, min(Eaaa), min(Baaa)]))
    Rx_ratio.set_ylim( 0,2)
    Rx_slope.set_xlim(0,0.4)
if 1:
    Rx_amp.set_ylabel('Amp') 
    Rx_slope.set_ylabel('Slope')
    Rx_ratio.set_ylabel('B/E')
    Rx_ratio.set_xlabel(fs.x_label())
    Rx_ratio.set_xscale(fs.x_scale)

    Rx_amp.set_yscale('log')

if 0:
    Rx_amp.set_ylim(  bump2(ext_amp, log=True))
    Rx_ratio.set_ylim(bump2(ext_ratio,log=False))
    Rx_slope.set_ylim(bump2(ext_slope,log=False))
if 1:
    Rx_amp.set_ylim(  bump2(ext_amp, log=True))
    Rx_ratio.set_ylim(0.0,2.0)
    Rx_slope.set_ylim(-10,0)
if 1:
    Rx_slope.set_xlim(bump2(ext_horz,log=False))
if len(figs) > 1:
    Rx_amp.set_xlim(bump2(ext_horz,log=False))
    Rx_ratio.set_xlim(bump2(ext_horz,log=False))
    Rx_slope.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    Rx_ratio.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    Rx_amp.set_xlabel(fs.x_label())
    Rx_slope.set_xlabel(fs.x_label())



for nfig,fig in enumerate(figs):
    fig.set_figwidth(13)
    outname = out_format%(car_outname,field,fig_names[nfig],plot_format)

    fig.savefig(outname)
    print outname
    plt.close(fig)
