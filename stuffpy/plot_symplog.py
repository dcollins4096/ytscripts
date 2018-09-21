def imshow_symlog(my_matrix, vmin, vmax, logthresh=5,this_plt=plt):
    img=this_plt.imshow( my_matrix ,
                vmin=float(vmin), vmax=float(vmax),
                norm=matplotlib.colors.SymLogNorm(10**-logthresh) )

    maxlog=int(np.ceil( np.log10(vmax) ))
    minlog=int(np.ceil( np.log10(-vmin) ))

    #generate logarithmic ticks 
    tick_locations=([-(10**x) for x in xrange(minlog,-logthresh-1,-1)]
                    +[0.0]
                    +[(10**x) for x in xrange(-logthresh,maxlog+1)] )

    cb=this_plt.colorbar(ticks=tick_locations)
    return img,cb

plt.clf()
#imshow_symlog(theta[1:-1,1:-1].v, -1.5, 1.5)
#plt.imshow(np.abs(theta[1:-1,1:-1].v),norm = matplotlib.colors.Normalize(0,1.5))
#plt.imshow(np.log10(MPhi[1:-1,1:-1].v))
#plt.imshow(MPhi[1:-1,1:-1].v, norm = matplotlib.colors.Normalize(1,1e3))
#plt.imshow(np.log10(frb_d[1:-1,1:-1].v))
#plt.imshow(np.log10(Phi[1:-1,1:-1].v))
#grad_n = grad_y*grad_y+grad_z*grad_z
#plt.imshow(np.log10(grad_n[1:-1,1:-1].v))
if 0:
    tmp = copy.copy(frb_d[1:-1,1:-1].v*1.4e22)
    tmp[tmp > 5e21] = 0
    plt.imshow(tmp)
    fname = 'p07_SigmaLT5.png'
    print fname
if 0:
    tmp = copy.copy(MPhi[1:-1,1:-1].v)
    tmp[tmp > 10] = 0
    plt.imshow(tmp)
    fname = 'p07_MPhiLT10.png'
    print fname
    plt.colorbar()
    plt.savefig(fname)

if 1:
#30*16x 25
    L = nar([(70.)/512,55./512,0])
    R = nar([(100.)/512   ,(55.+25.)/512,1])
