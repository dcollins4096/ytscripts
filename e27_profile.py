
basedir = '/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/'
b01 = 'b01_ct_turb_timing'
b02 = 'b02_ded_turb_timing'

fname1 = "%s/%s/%s"%(basedir,b01,'gprof_stock_hlld_hack.txt')
fname2 = "%s/%s/%s"%(basedir,b02,'gprof_stock_hacked.txt')

fptr1 = open(fname1)
times1 = []
calls1 = []
for line in fptr1:
    spl = line.split('!')
    T = float(spl[0])
    F = spl[1]
    times1.append(T)
    calls1.append(F)
fptr2 = open(fname2)
times2 = []
calls2 = []
for line in fptr2:
    spl = line.split('!')
    T = float(spl[0])
    F = spl[1]
    times2.append(T)
    calls2.append(F)

plt.clf()
if 0:
    plt.plot(np.cumsum(times1),'r',label='ct', marker="+")
    plt.plot(np.cumsum(times2),'g',label='ded', marker="+")
    outname = 'times_cuml.pdf'
if 0:
    plt.plot(np.cumsum(times1[::-1])[::-1],'r',label='ct')
    plt.plot(np.cumsum(times2[::-1])[::-1],'g',label='ded')
    outname = 'times_cuml_rev.pdf'
if 1:
    plt.plot(times1,'r',label='ct',marker="+")
    plt.plot(times2,'g',label='ded',marker="+")
    outname = 'times.pdf'

plt.xlabel('Call number')
plt.ylabel('time per call')
plt.legend(loc=1)
plt.savefig(outname)
print outname
