from go import *

def gaa(x,sigma=None):
    original_curve= np.exp(-x**2/(2*sigma**2))
    original_curve/=original_curve.sum()
    return original_curve
sigma=1.
x = np.arange(-20*sigma, 20*sigma, 1/256)
dx=x[1]-x[0]
curve_g= gaa(x,sigma=sigma)
if 1:
    curve_p= (1+(x/(sigma/3))**2)**(-0.5)
    curve_p*= curve_g.max()/curve_p.max()
if 0:
    a=np.pi/2
    curve_p = 1/(np.exp(a*x)+np.exp(-a*x))
    curve_p*= curve_g.max()/curve_p.max()
if 0:
    a=np.pi/2
    curve_p = np.exp(x)/(1+np.exp(x))**2
    curve_p*= curve_g.max()/curve_p.max()
sigma_p= ((x**2*curve_p*dx).sum()/(curve_p*dx).sum() )**0.5
#print("var %.2f %.2f"%(sigma_p, np.pi/2/a))
fig,ax=plt.subplots(1,3,figsize=(12,8))
ax[0].plot(x,curve_p,c='g')
for nc,original_curve in enumerate([curve_g,curve_p]):

    orig_sigma= ((x**2*original_curve*dx).sum()/(original_curve*dx).sum() )**0.5
    print(orig_sigma)
    ax[nc].plot(x,original_curve,c='k')
    widths=[]
    kerns=np.linspace(0.5,2,20)
    #kerns=nar([0.5])
    for kern in kerns:
        gaussian = gaa(x,sigma=kern)
        result = np.convolve(original_curve, gaussian, mode="same")
        wid = (x**2*result*dx).sum()/(result*dx).sum()
        widths.append( wid)

        ax[nc].plot(x,result,c='b')

    ax[2].scatter(kerns, np.sqrt(widths-kerns**2))
    #ax[2].plot(kerns, 1-np.sqrt(widths-kerns**2)/orig_sigma)
#ax[2].set_ylim(sigma_p-0.1,sigma_p+0.1)
fig.savefig(os.environ['HOME']+'/PigPen/p44_widths.png')
