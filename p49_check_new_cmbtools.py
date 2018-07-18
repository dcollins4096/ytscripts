import astropy.io.fits as pyfits


basedir = '/Users/dcollins/RESEARCH2/Paper49_EBQU/2018-05-21-EBQU-verify/NASA_test_cmbtools/ab19_new'
plt.clf()

stuff={}
stuff['orig_stampede']={'dir':'FRB_orig_stampede','temp':'DD%04d_%s%s_n0-0001_p-1.fits'}
stuff['orig_pfe']={'dir':'FRB_orig_pfe','temp':'DD%04d_%s%s.fits'}
stuff['new_pfe']={'dir':'FRB_new_pfe','temp':'DD%04d_%s%s.fits'}
page = 'stuff.html'
hptr = open(page,'w')
frame=275
sim='ab19'
hptr.write('<table>')
for suite in ['new_pfe','orig_pfe','orig_stampede', 'density']:
    hptr.write('<tr><th> %s </th>\n'%suite)
    for quan in  'EBQU':
        hptr.write('<td style="background-color:yellow;"> %s </td>'%quan)
        for ax in 'xyz':
            if suite not in ['density']:
                fname = "%s/%s/%s"%(basedir,stuff[suite]['dir'],stuff[suite]['temp']%(frame,quan,ax))
                plt.clf()
                data=np.array(pyfits.open(fname)[0].data,dtype=np.double)
                plot = plt.imshow(data,interpolation='nearest',origin='lower')
                plt.colorbar(plot)
                plt.title("%s %s %s"%(suite,quan,ax))
                outname = 'p19_check_tool_%s_%s_%04d_%s_%s.png'%(sim,suite,frame,quan,ax)
                hptr.write('<td>')
                hptr.write('<img src=%s>'%outname)
                hptr.write('</td>\n')
                plt.savefig(outname)
                print(outname)
            if suite in ['density']:
                fname = "%s/%s/%s"%(basedir,stuff['new_pfe']['dir'],stuff['new_pfe']['temp']%(frame,'density','_'+ax))

                data=np.array(pyfits.open(fname)[0].data,dtype=np.double)
                plt.clf()
                plot = plt.imshow(data,interpolation='nearest',origin='lower')
                plt.title("%s %s %s"%('new_pfe','density',ax))
                plt.colorbar(plot)
                hptr.write('<td>')
                outname = 'p19_check_tool_%s_%s_%04d_%s_%s.png'%(sim,'new_pfe',frame,'density',ax)
                hptr.write('<img src=%s>'%outname)
                hptr.write('</td>\n')
                plt.savefig(outname)
                print(outname)
    hptr.write('</tr>\n')
hptr.write('</table>\n')
hptr.close()
