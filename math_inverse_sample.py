
#mooch from here:http://matlabtricks.com/post-44/generate-random-numbers-with-a-given-distribution
from scipy.interpolate import interp1d
from scipy import interpolate
plt.clf()
P = np.random.uniform(30)

prob = nar([0, 1, 10, 2, 0, 0, 4, 5, 3, 1, 0])
x = np.arange(len(prob))
x_fine = np.linspace(0,x.max(),100)
#pdf = interp1d( x, prob, kind='cubic')
spline_tool = interpolate.splrep(x,prob)
pdf = interpolate.splev(x_fine, spline_tool)
pdf=pdf/pdf.sum()
neg = pdf<0
pdf[neg]=0
#pdf = interpolate.( x, prob, kind='cubic')
#plt.plot(x_fine,pdf(x_fine),c='r')
plt.plot(x_fine, pdf, c='g')
plt.scatter(x,prob,c='b')
plt.savefig('math_probest.pdf')

#make a cdf
cdf = np.cumsum(pdf)/cdf.max()
val,ind = np.unique(cdf,return_index=True)
plt.clf()
plt.scatter(x_fine[ind],val,c='r')
plt.savefig('math_probest_cdf2.pdf')
dumb_plt(plt,x_fine,cdf,'x','cdf','math_probtest_cdf.pdf')

#now a random thing
tool=interp1d(val,x_fine[ind])
should_work = tool(np.random.uniform(size=1200))
plt.clf()
#oot=plt.hist(should_work,histtype='step',normed=True,bins=100)
oot=np.histogram(should_work,bins=100)
bincenters=0.5*(oot[1][1:]+oot[1][:-1])
bins = oot[0]*1.0/oot[0].sum()
plt.plot(bincenters,bins,c='r')
plt.savefig('math_probest_hist.pdf')
dumb_plt(plt,x_fine,pdf,'x','pdf','math_probtest_with_hist.pdf')
