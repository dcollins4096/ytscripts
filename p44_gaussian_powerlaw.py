if 'ef' not in dir():
    execfile('go')
#from p44_random_tools import *
ef('p44_random_tools.py')
#plt.clf()
#r1 = np.random.normal(0,1,size=1000)
#plt.hist(r1,histtype='step',bins=30)
#plt.savefig('p44_power2.pdf')
"""
double Gaussian(double cs)
{

  double mean = 0;
  double stdev = cs;
  double u1 = rand();
  u1 = u1/RAND_MAX;
  double u2 = rand();
  u2 = u2/RAND_MAX;
  double x = mean + stdev*sqrt(-2*log(u1))*cos(2*M_PI*u2);

  return x;
}
     k_wave = sqrt(k2)*2.0*pi;
     Ak0 = Gaussian(pow(k_wave, -ind*0.5-1));
     Ax Ay Az =  random unit vector
          Az /= RAND_MAX;
          AA = sqrt(Ax*Ax+Ay*Ay+Az*Az);
          Ax /= AA;
     vel[0][igrid] += Ak0*(-ky*Az*sin(kdotx+phiz) + kz*Ay*sin(kdotx+phiy));
     vel[1][igrid] += Ak0*(-kz*Ax*sin(kdotx+phix) + kx*Az*sin(kdotx+phiz));
     vel[2][igrid] += Ak0*(-kx*Ay*sin(kdotx+phiy) + ky*Ax*sin(kdotx+phix));

"""

V = []
slopes=[]
powers=[]

Nx = 3000
dk = 2*np.pi/Nx
k = np.arange(0,2*np.pi-0.2*dk,dk)
slope = 0.2
n_runs = 1
amps=[]
ap = amp_phase(k,slope)
for n in range(n_runs):
    amp = ap(a='D',p='U')
    amps.append(ap.amp*ap.phases)
    Vhat = symmetric(amps[-1])
    V.append(np.fft.ifft(Vhat))

plt.clf()
plt.hist(ap.amp,histtype='step')
plt.ylabel='N'
plt.xlabel='K amplitude'
plt.savefig('p44_ap_amp.pdf')

fig, ((ax_power, ax_slope), (ax_pdf, ax_pdf_k)) = plt.subplots(2, 2)
fig.tight_layout()

format1 = "%10s N %f x0 %f sigma %f"
ax_pdf.set_title('x-space pdf')
all_v=nar(V).flatten().real
val,bins=np.histogram(all_v, bins=50)
b = 0.5*(bins[1:]+bins[:-1])
val = 1.0*val/val.sum()
ax_pdf.plot(b,val)
v_fit = gauss_fit(b,val)
ax_pdf.plot(b,gauss_me(b, v_fit['fit_norm'],v_fit['fit_center'],v_fit['fit_width']))
print format1%("V fit", v_fit['fit_norm'],v_fit['fit_center'],v_fit['fit_width'])
ax_pdf.set_xlabel('v')
ax_pdf.text(-0.1,0,r"$v_0=%0.2e \sigma = %0.2e$"%(v_fit['fit_center'], v_fit['fit_width']))


ax_pdf_k.set_title('k-space pdf')
all_a=nar(amps).flatten().real
val,bins=np.histogram(all_a, bins=50)
b = 0.5*(bins[1:]+bins[:-1])
val = 1.0*val/val.sum()
ax_pdf_k.plot(b,val,'b-')
ar_fit = gauss_fit(b,val)
print format1%("A real fit", ar_fit['fit_norm'],ar_fit['fit_center'],ar_fit['fit_width'])
ax_pdf_k.plot(b,gauss_me(b, ar_fit['fit_norm'],ar_fit['fit_center'],ar_fit['fit_width']), 'r-')

all_a=nar(amps).flatten().imag
val,bins=np.histogram(all_a, bins=50)
b = 0.5*(bins[1:]+bins[:-1])
val = 1.0*val/val.sum()
#ax_pdf_k.plot(b,val,'b--')
ai_fit = gauss_fit(b,val)
print format1%("A imag fit", ai_fit['fit_norm'],ai_fit['fit_center'],ai_fit['fit_width'])
ax_pdf_k.plot(b,gauss_me(b, ai_fit['fit_norm'],ai_fit['fit_center'],ai_fit['fit_width']),'r--')
ax_pdf_k.set_xlabel('v')
ax_pdf_k.text(-0.1,0,r"$v_0=%0.2e \sigma = %0.2e$"%(ai_fit['fit_center'], ai_fit['fit_width']))


#ax_pdf.plot(k,np.exp(-0.5**(k-k/2)**2))
avg_power=0
avg_slope=0
for n in range(n_runs):
    slope_i,kmin=check_powerlaw(V[n])
    slopes.append(slope_i)
    power_i = power_spectrum(V[n])
    powers.append(power_i)
    avg_power = power_i + avg_power
    avg_slope = slope_i + avg_slope
    ax_slope.plot(k[1:slope_i.size],slope_i[1:],c=[0.5]*3)
    ax_power.plot(k[1:power_i.size],power_i[1:],c=[0.5]*3)

dumb_plt(ax_slope,k[kmin+1:Nx/2],nar(avg_slope)/len(V), 'k','<alpha>','P44_gaussian_powerlaw.pdf',c='k')
dumb_plt(ax_power,k[1:Nx/2],nar(avg_power[1:Nx/2])/len(V), 'k','<Pf>','P44_gaussian_powerlaw.pdf',scale=('linear','log'),c='k')
dumb_plt(ax_power,k[1:Nx/2],k[1:Nx/2]**(-slope), 'k','<Pf>','P44_gaussian_powerlaw.pdf',scale=('linear','log'),c='r')
plt.close(fig)
