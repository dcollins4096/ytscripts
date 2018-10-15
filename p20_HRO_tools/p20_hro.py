import pyfits

hsv_copy = copy.copy(matplotlib.cm.hsv)

basedir = '/Users/dcollins/RESEARCH2/Paper20_Filaments'

B1 = 'B2Hpix_masked.fits'
f1 = pyfits.open('%s/%s'%(basedir,B1))
f1d = f1[0].data
masked_array = np.ma.array(f1d, mask=np.isnan(f1d))
hsv_copy.set_bad('k',1)
plt.clf()
plt.imshow(masked_array,origin='lower',interpolation='nearest',cmap=hsv_copy)
plt.savefig('test1.png')
