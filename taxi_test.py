
dirlist={}
base = '/scratch1/dcollins/Paper08/'
dirlist['b02_5'] = 'B02/512'
dirlist['b02_2'] = 'B02/256'
dirlist['b20_2'] = 'B20/256'
dirlist['b2_2'] = 'B2/256'
dirlist['b20_5'] = 'B20/512'
dirlist['b2_5'] = 'B2/512'
framelist = [10,50]
all_x=['Bstretching', 'Bcompression', 'Badvection', 'Ltension', 'Lpressure']

labels={}
Fcol= {'Bstretching':'r', 'Bcompression':'g', 'Badvection':'b', 'Ltension':'c', 'Lpressure':'m'}
Line={}
for sim in ['b2_2']:
    for frame in framelist:
        plt.clf()
        ds = yt.load('%s/%s/RS%04d/restart%04d'%(base,dirlist[sim],frame,frame))
        if 1:
            """for tsting"""
            if 'peak' not in dir():
                value, peak = ds.find_max('density')
            width = 0.05
            data = ds.h.sphere(peak,width)
        if 0:
            data = ds.all_data()
        for field_base in all_x:
            flavor = "_f"
            field = field_base+flavor
        
            prof = yt.create_profile(data,'density',field,weight_field='cell_volume',logs={field:False})
            field_data = prof.field_data
            prof.field_data['x_bin']=prof.x_bins
            fff = prof.field_data[('gas',field)]
            plt.plot(0.5*(field_data['x_bin'][1:]+field_data['x_bin'][:-1]), fff,linestyle = Line.get(field,'-'),
                     label=labels.get(field,field), c=Fcol.get(field_base,'r'))
            plt.xscale('log'); plt.yscale('symlog')
            plt.xlabel('density'); plt.ylabel(field)
        outname = '%s_n%04d_%s_%s.png'%(sim,frame,'stuff',flavor)
        plt.legend(loc=0)
        plt.savefig(outname)
        print outname



