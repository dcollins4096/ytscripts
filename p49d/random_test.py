
z,y,x = np.mgrid[-1:1:10j,-1:1:10j,-1:1:10j]
if 'ra' not in dir():
    ra = np.random.random(z.shape)
r2 = x*x+y*y+z*z

ok = np.zeros_like(r2)
ok[r2<0.5] = 1
ok += 0.1*ra
mok=np.mean(ok)
mok2 = np.sum(ok)/ok.size
sok = np.std(ok)
sok2 = np.sqrt(np.sum( (ok-mok)**2)/ok.size)
sok3 = np.sqrt( np.sum(ok**2 - 2*ok*mok + mok**2) /ok.size) 
sok4 = np.sqrt( np.sum(ok**2 ) / ok.size + np.sum( (- 2*ok*mok + mok**2) /ok.size) )
sok5 = np.sqrt( np.sum(ok**2 ) / ok.size - mok**2 ) # np.sum( (- 2*ok*mok + mok**2) /ok.size) )
print(sok5-sok)
