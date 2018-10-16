ef('p44_random_tools.py')

fig4, ((ax_k,ax_x),(ax_vk,ax_vx),(ax_pk,ax_px)) = plt.subplots(3, 2)
#ax_k.plot(k,np.sin(k))
#ax_k.plot(k,np.abs(vhat))
if 1:
    form2 = "%10s %f"
    x_width=[]
    k_width=[]
    rat=[]
    nx = []

    #nx_list = range(1000,8000,100)
    fourposter = True
    gaussfit = True
    #nx_list=[10000]
    nx_list=[1000]
    #nx_list = [7000] # 
    #sigma_list = np.arange(0.5,10,0.5)
    sigma_list = [3]
    fname = 'p44_4panel.pdf'
    for Nx in nx_list: #Nx = 3000
        for sigma in sigma_list:
            dk = 2*np.pi/Nx
            k = np.arange(0,2*np.pi-0.2*dk,dk)
            if 0:
                ap = amp_phase(k,0.2)
                ap.k=k
                ap(a='B',p='U') #,sigma=sigma )
            Ak = ap.amp*ap.phases
            #vhat = Ak #symmetric(Ak)
            vhat = symmetric(Ak)
            if 1:
                #print "var vhat", (np.abs(vhat)**2).sum()/Nx
                #print "var v   ", np.var(v)
                vhat[0] = 0.0 # 3.0*Nx # (sigma*1.0)*Nx
            v = np.fft.ifft(vhat)  *np.sqrt(Nx)
            if np.abs( np.abs(v.imag).sum()/np.abs(v.real).sum() ) > 1e-10:
                print "REAL ERROR"
            v=v.real
            if 1:
                print "var vhat", (np.abs(vhat)**2).sum()/Nx
                vbar= np.sum(v)/Nx
                print "var v   ", (np.abs(v-vbar)**2).sum()/Nx

            if fourposter: 
                scale=('log','log')
                sl = slice(1,None)
                dumb_plt(ax_px, k[:Nx/2][sl], (np.abs(np.fft.fft(v)[:Nx/2])**2)[sl],'','Px',fname,scale=scale)
                dumb_plt(ax_pk, k[sl], (np.abs(vhat)**2)[sl],'','Pk',fname,scale=scale)

            if fourposter: dumb_plt(ax_k,k,vhat, '','vhat','p44_4panel.pdf')
            #ax_k.set_ylim(0,1.1)
            if fourposter: dumb_plt(ax_x,k,v,'','v',fname)
            #print form2%("imag", np.sum( np.abs( v.imag)))

            v_val,v_bins=np.histogram(v, bins=50)
            v_b = 0.5*(v_bins[1:]+v_bins[:-1])
            #v_val = 1.0*v_val/v_val.sum()
            if gaussfit: 
                    v_fit = gauss_fit(v_b,v_val)
                    x_width.append(v_fit['fit_width'])
            if fourposter: 
                if gaussfit:
                    print form2%("x sig",v_fit['fit_width']), v.shape, "target",sigma
                    if fourposter: ax_vx.plot(v_b,gauss_me(v_b, v_fit['fit_norm'],v_fit['fit_center'],v_fit['fit_width']),'r--')
                dumb_plt(ax_vx,v_b,v_val,'','V(x)',fname,c='r')

            k_val,k_bins=np.histogram(vhat.real, bins=50)
            k_b = 0.5*(k_bins[1:]+k_bins[:-1])
            k_val = 1.0*k_val/k_val.sum()
            if gaussfit:
                k_fit = gauss_fit(k_b,k_val)
                k_width.append(k_fit['fit_width'])
                g_should = gauss_me(k_b, 1,0,sigma**2)
                g_should/=g_should.sum()
                print form2%("k sig",k_fit['fit_width']), k.shape
                if fourposter:
                    ax_vk.plot(k_b,gauss_me(k_b, k_fit['fit_norm'],k_fit['fit_center'],k_fit['fit_width']),c='g')
                    ax_vk.plot(k_b,g_should,c='k')
            if fourposter: dumb_plt(ax_vk,k_b,k_val,'','V(k)',fname,c='b')

            if gaussfit:
                rat.append(k_fit['fit_width']/v_fit['fit_width'])
            nx.append(Nx)

fig4.savefig('p44_4panel.pdf')
plt.close(fig4)

plt.savefig('p44_normal_in_k.pdf')
if 1:
    nx_list = nar(nx_list)
    rm = rainbow_map(max(nx_list))
    fig2, ax2 = plt.subplots(1, 1)
    x_width=nar(x_width)
    k_width=nar(k_width)
    rat=nar(rat)
    nx=nar(nx)
    #dumb_plt(ax2,k_width/np.sqrt(nx),x_width,'Sk/sqrt(n)','Sx','p44_fft_gauss_2.pdf')
    #dumb_plt(ax2,k_width,x_width,'Sk','Sx','p44_fft_gauss_2.pdf')
    dumb_plt(ax2,k_width,x_width,'k_width','x_width','p44_fft_gauss_4.pdf',scatter=True,c=rm(nx))
    #ax2.scatter(k_width,k_width/x_width,c=rm(nx_list))
    #fig2.savefig('p44_fft_gauss_4.pdf')
    plt.close(fig2)
if 0:
    the_x,the_y = nx,(x_width/k_width)**2
    x0,x1=the_x.min(),the_x.max()
    fit = np.polyfit(the_x, the_y, 1)
    m = fit[0]
    b = fit[1]
    def the_f(x):
        return m*x+b
    ax2.scatter(the_x,the_y,c=rm(nx))
    ax2.plot([x0,x1],[the_f(x0),the_f(x1)])
