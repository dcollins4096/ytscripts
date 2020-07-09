
from go import *
import turb_quan
reload(turb_quan)
proto_list_very_high_field=['za01','zb01','zc01','zd01','ze01']
proto_list_Mach10=['za05','zb05','zc05','zd05']#,'ze01']
proto_list_Mach1=['za04','zb04','zc04','zd04']#,'ze01']
proto_list_256=['zc01','zc04','zc05']
proto_list_512=['zd01','zd04','zd05']
box_list_standard = [1,2,  4, 8, 16, 32]
box_list_256 = [4,4,4]
box_list_512 = [8,8,8]
SERIES = 'zX05'

axlist = 'x'

box_list = box_list_standard
car_list = proto_list_Mach10
car_list = ['ca02']
box_list=[1]
#car_list=['zc01']; box_list=[4]
#box_list = [1, 2, 2, 2, 1]
#car_list=['za01']#,'zc01']
#car_list=['zc01']#,'zd01','ze01']
#car_list=['za02','zb02','zc02','zd02']#,'ze02']
frames_t01 = {'ze01':127, 'zd01':11,'za01':11,'zb01':11,'zc01':11}
frames_256 = {'zc01':31,'zc04':101,'zc05':50}
frames_512 = {'zd01':16,'zd04':101,'zd05':16}
car_array=[]
for name in car_list:
    car_array.append(taxi.load(name))
    car_array[-1].frames=list(range(10))+list(range(10,150,10))

    #car_array[-1].frames=[10,13,16,31] #[100]
    #car_array[-1].frames=[7,9,16] #[10,13,16,31] #[100]
    #car_array[-1].frames=[frames_t01[name]]
    #car_array[-1].frames=[ frames_256[name]]
    #car_array[-1].frames=[ frames_512[name]]


if 1:
    for nc,car in enumerate(car_list):
        car_array[nc].box_size = box_list[nc]

flt = taxi.fleet(car_array)
plt.close('all')
def make_eb(carname,frames=None):
    qb = turb_quan.quan_box(car=taxi.load(carname))
    qb.EBall(frames=frames)

#for carname in car_list:
#    make_eb(carname,frames=[7,9,16])
#make_eb('za03',frames=[101])
#make_eb('zc02',frames=[101])
#make_eb('zd02',frames=[101])

if 0:
    for nax,ax in enumerate(axlist):
        for nf,frame in enumerate(car_array[0].frames):
            d={};h={};q={};u={};e={};b={}
            for car in flt.taxi_list:
                set_output = "TD%04d"%nf #for when the DD outputs are irregular.
                set_output = "DD%04d"%frame
                xd='DD'
                if car.name in ['zd01','ze01']:
                    xd='XD'
                density= "%s/frbs/%s%04d_density_%s.fits"%(car.directory,xd,frame,ax)
                field= "%s/frbs/%s%04d_magnetic_field_strength_%s.fits"%(car.directory,xd,frame,ax)
                Q= "%s/frbs/%s%04d_Q%s.fits"%(car.directory,xd,frame,ax)
                U= "%s/frbs/%s%04d_U%s.fits"%(car.directory,xd,frame,ax)
                E= "%s/frbs/%s%04d_E%s.fits"%(car.directory,xd,frame,ax)
                B= "%s/frbs/%s%04d_B%s.fits"%(car.directory,xd,frame,ax)
                for fname in [Q,U,E,B]:
                    if not os.path.exists(fname):
                        print("shoot.  Missing %s"%fname)
                        continue
                #if not os.path.exists(field):
                #    print("shoot.  Missing %s"%field)
                d[car.name]=di=pyfits.open(density,dtype=np.double)[0].data
                h[car.name]=hi=pyfits.open(field,dtype=np.double)[0].data
                q[car.name]=qi=pyfits.open(Q,dtype=np.double)[0].data
                u[car.name]=ui=pyfits.open(U,dtype=np.double)[0].data
                e[car.name]=ei=pyfits.open(E,dtype=np.double)[0].data
                b[car.name]=bi=pyfits.open(B,dtype=np.double)[0].data
                if 1:
                    outname = "p49six_%s_%s_%s.png"%(car.outname,set_output,ax)
                    turb_quan.plotter2([d[car.name], h[car.name],q[car.name],
                                        u[car.name],e[car.name],b[car.name]],
                                           outname,
                                           labs=['d','H','q','u','e','b'],
                                           norm='ind')
                    if 0:
                        outname = 'p49_qhist_%s_%s.png'%(car.outname,frame,set_output)
                        plt.hist(np.log(np.abs(q.flatten())),histtype='step',color='rgb'[nax])
                if 1:
                    outname_base = "%s__%04d_frb_%s_%s.png"
                    for name, array in  [ ['Q',qi],['U',ui],['E',ei],['B',bi]]:
                        outname = outname_base%( car.outname, frame, name, ax)
                        turb_quan.plotter2([array], outname, labs=[name], norm='ind',npx=1)

                
            if 0:
                outname = "p49six_density_%s_%s_%s.png"%(SERIES,set_output,ax)
                turb_quan.plotter2( [d[name] for name in car_list],
                                    outname,share=False,
                                       labs=car_list,
                                       norm='symlog')
            if 0:
                outname = "p49six_hfield_%s_%s_%s.png"%(SERIES,set_output,ax)
                turb_quan.plotter2( [h[name] for name in car_list],
                                    outname,share=False,
                                       labs=car_list,
                                       norm='symlog')
                outname = "p49six_Q_%s_%s_%s.png"%(SERIES,set_output,ax)
                turb_quan.plotter2( [q[name] for name in car_list],
                                    outname,share=False,
                                       labs=car_list,
                                       norm='symlog')
                outname = "p49six_U_%s_%s_%s.png"%(SERIES,set_output,ax)
                turb_quan.plotter2( [u[name] for name in car_list],
                                    outname,share=False,
                                       labs=car_list,
                                       norm='symlog')
                outname = "p49six_E_%s_%s_%s.png"%(SERIES,set_output,ax)
                turb_quan.plotter2( [e[name] for name in car_list],
                                    outname,share=False,
                                       labs=car_list,
                                       norm='symlog')
                outname = "p49six_B_%s_%s_%s.png"%(SERIES,set_output,ax)
                turb_quan.plotter2( [b[name] for name in car_list],
                                    outname,share=False,
                                       labs=car_list,
                                       norm='symlog')

#import cmbtools
import p49_QU2EB
reload(p49_QU2EB)
def Cl2(Q,U,Density,car=None):
    stuff=p49_QU2EB.EBfromQU(Q,U,BoxSize=car.box_size,return_quharm=True)
    Eharm=stuff['Eh']
    Bharm=stuff['Bh']
    Deltal = stuff['Deltal']
    Delta=stuff['Delta']
    N=stuff['N']
    lmax = Deltal[0]*N[0]
    lbins = np.linspace(0,lmax,64*car.box_size/2)
    lcent = lbins[:-1] + np.diff(lbins)/2.
    ClEE = cmbtools.harm2cl(Eharm,Deltal,lbins)
    ClBB = cmbtools.harm2cl(Bharm,Deltal,lbins)
    stuff.update( {'ClEE':ClEE,'ClBB':ClBB,'lbins':lbins,'lcent':lcent})
    return stuff
    #ClTE = cmbtools.harm2clcross_samegrid(Eharm,Nharm,Deltal,lbins)
    #return {'clee':ClEE,'clbb':ClBB,'ell':lbins,'delta':Delta,'deltal':Deltal,
    #        'Eh':Eharm,'Bh':Bharm,'Qh':Qharm,'Uh':Uharm,'Q':Q,'U':U,'xsize':xsize}

if 1:
    fig_cl, ax_cl = plt.subplots(1,1)
    wtfmate=[]
    #all sims on one plot
    fixed_frames=car_array[0].return_frames()
    stuff={}
    for car in flt.taxi_list:
        for nf,frame in enumerate(car.return_frames()):
            set_output = "DD%04d"%frame
            d,h,q,u,e,b= [],[],[],[],[],[]
            for nax,ax in enumerate(axlist):
                xd = 'DD'
                if car.name in ['zd01','ze01']:
                    xd='XD'
            #set_output = "TD%04d"%nf #for when the DD outputs are irregular.
                density= "%s/frbs/%s%04d_density_%s.fits"%(car.directory,xd,frame,ax)
                field= "%s/frbs/%s%04d_magnetic_field_strength_%s.fits"%(car.directory,xd,frame,ax)
                Q= "%s/frbs/%s%04d_Q%s.fits"%(car.directory,xd,frame,ax)
                U= "%s/frbs/%s%04d_U%s.fits"%(car.directory,xd,frame,ax)
                E= "%s/frbs/%s%04d_E%s.fits"%(car.directory,xd,frame,ax)
                B= "%s/frbs/%s%04d_B%s.fits"%(car.directory,xd,frame,ax)
                if not os.path.exists(E):
                    print("shoot.  Missing %s"%E)
                    continue
                #if not os.path.exists(field):
                #    print("shoot.  Missing %s"%field)
                d.append(np.ascontiguousarray(pyfits.open(density)[0].data,dtype=np.double))
                #h.append(np.ascontiguousarray(pyfits.open(field,dtype=np.double)[0].data))
                q.append(np.ascontiguousarray(pyfits.open(Q)[0].data,dtype=np.double))
                u.append(np.ascontiguousarray(pyfits.open(U)[0].data,dtype=np.double))
                e.append(np.ascontiguousarray(pyfits.open(E)[0].data,dtype=np.double))
                #b.append(np.ascontiguousarray(pyfits.open(B,dtype=np.double)[0].data))
                if 1:
                    import cmbtools
                    #x = cmbtools.map2harm(Q,np.ones(2))
                    ts=Cl2(q[-1],u[-1],d[-1],car)
                    stuff[car.name]=ts
                    labial = None
                    labial = "%s %d"%(car.name ,frame)
                    linestyle = {10:'-',13:'--',16:':',31:'-.',101:'-'}.get(frame,'-')
                    cdict = {'za01':'r','zb01':'g','zc01':'b','zd01':'k','ze01':'c'}
                    cdict.update({'za05':'r','zb05':'g','zc05':'b','zd05':'k','ze05':'c'})
                    cdict.update({'za04':'r','zb04':'g','zc04':'b','zd04':'k','ze04':'c'})
                    #cdict['zc01']='r';cdict['zc04']='g';cdict['zc05']='b'
                    #cdict['zd01']='r';cdict['zd04']='g';cdict['zd05']='b'
                    c=cdict[car.name]
                    print('smooth')
                    kvals = ts['lcent']
                    comp = 0#2.35
                    #ax_cl.plot(kvals,ts['ClBB']*kvals**comp, label=labial ,linestyle=linestyle,c=c)
                    #ax_cl.plot(kvals,ts['ClEE']*kvals**comp, linestyle='-',c=c,label=labial)
                    #ax_cl.plot(kvals,ts['ClBB']*kvals**comp, linestyle='--',c=c)
                    #ax_cl.plot(kvals,ts['ClBB']/ts['ClEE'], linestyle='--',c=c)
                    ax_cl.plot(kvals,ts['ClBB']/ts['ClEE'], linestyle='--',label=labial)
                    wtfmate.append(car.name)
    #powerline(ax_cl,ts['lcent'][3], ts['lcent'][-1], 100*ts['ClEE'][0],-2.5)
    ax_cl.legend(loc=0)
    ax_cl.set_xscale('log'); ax_cl.set_yscale('log')
    ylabel = r'$C_\ell^{BB}$'
    ylabel = r'$C_\ell^{EE}$'
    if np.abs(comp) > 1e-5:
        ylabel += r'$k^{%.2f}$'%comp
    ylabel = r'$C_\ell^{BB}/C_\ell^{EE}$'
    ax_cl.set_xlabel('k'); ax_cl.set_ylabel(ylabel)
    outname = 'p49_multispectra_EEBB__%s_%s_16.png'%(ax,SERIES)

    fig_cl.savefig(outname)
    print(outname)
    if 0:


        if 0:
            outname = "p49bysim_%s_%s_%s_%s.png"%(car.outname,set_output,'density',ax)
            turb_quan.plotter2(d,
                               outname,
                               labs=[car.outname for car in flt.taxi_list],
                               norm='ind',share=False)
        if 0:
            outname = "p49bysim_%s_%s_%s_%s.png"%(car.outname,set_output,'density',ax)
            turb_quan.plotter2(d,
                               outname,
                               labs=[car.outname for car in flt.taxi_list],
                               norm='ind',share=False)
        if 0:
            outname = 'p49_qhist_%s_%s.png'%(car.outname,frame,set_output)
            plt.hist(np.log(np.abs(q.flatten())),histtype='step',color='rgb'[nax])
            #print(outname)
            #plt.savefig(outname)
                    

if 0:
    fig_multi,ax_multi = plt.subplots(1,1)
    fig_ind,ax_ind = plt.subplots(1,1)
    thisax = ax_multi
    ell_list=[]

    for car in flt.taxi_list:
        xd = 'DD'
        if car.name in ['zd01','ze01']:
            xd='XD'
        qb = turb_quan.quan_box(car)
        for frame in car.frames:
            if 0: #check for files.
                for nax,ax in enumerate(axlist):
                    if 1: 
                        testfile= "%s/frbs/%s%04d_Q%s.fits"%(car.directory,xd,frame,ax)
                        if os.path.exists(testfile):
                            print("OK MAN %s"%testfile)
                        else:
                            print("Missing %s"%testfile)
                        continue
                continue
            #qb.EBall([frame])
            #qb.EBslopes([frame])
            
            for nax,ax in enumerate(axlist):

                clfile = "%s/frbs/%s%04d_Cl%s.dat"%(car.directory,xd,frame,ax)
                b = np.loadtxt(clfile)
                ell = b[:,0]
                ell_list.append(ell)
                EE = b[:,1]
                BB = b[:,2]
                ok = BB>0
                c='rgb'[nax]
                do=1
                if do==0:
                    thisax.plot( ell[ok], EE[ok]/BB[ok],c=c)
                    outname = "%s_DD%04d_ratio_take1.png"%(car.outname,frame)
                    thisax.set_ylim(1,5)
                    thisax.set_yscale('linear')
                if do==1:
                    thisax.plot( ell, EE, c=c,linestyle='-')#ell**3.8
                    thisax.plot( ell, BB, c=c,linestyle=':')#ell**3.8
                    outname = "%s_DD%04d_compensated_take1.png"%(car.outname,frame)
                    ylabel = r'$C_{xx}$'# \ell^{3.8}$'
                    #plt.ylim(5e2,1e4)
                    thisax.set_yscale('log')
            print(outname)
            thisax.set_xlabel(r'$\ell$')
            thisax.set_xscale('log')
            thisax.set_ylabel(ylabel)
            #fig_ind.savefig(outname)
    fig_multi.savefig('multiplots.png')

            #qb.QUEB(frame)
            #qb.EBslopes(frame)
