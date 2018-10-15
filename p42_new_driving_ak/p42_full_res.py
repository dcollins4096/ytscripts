


simname = 'B20'
frame = 69
output_dir_name = 'UnigridExtractions'
output_prefix = 'TwoLevel_Density'
basedir1 = '/scratch1/dcollins/Paper08/%s/512'%simname
out_dir = "%s/%s"%(basedir1,output_dir_name)
out_name = "%s/%s_%s_512_%04d"%(out_dir,output_prefix,simname,frame)
if 'fptr' in dir():
    del fptr
if 'den' not in dir():
    fptr = h5py.File(out_name,'r')
field = 'Density'
axis = 2
image_name = 'Projection_%s_%s_%s_%s_t%04d'%('xyz'[axis],output_prefix,simname, field, frame)
try:
    if 'den' not in dir():
        den3 = fptr[field][:]
        den = np.sum( fptr[field], axis=axis)
    plt.clf()
    if 1:
        fig = plt.figure(figsize=[20.48]*2, dpi=100)
        fig.figimage(np.log10(den))
        fig.savefig(image_name)
        print image_name
    if 0:
        plt.imshow(np.log10(den))
        plt.savefig(image_name)
        print image_name
except:
    raise
finally:
    if 'fptr' in dir():
        fptr.close()
