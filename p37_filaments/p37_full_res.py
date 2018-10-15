
plt.clf()
fig = plt.figure(figsize=[81.92]*2, dpi=100)
fig.figimage(np.log10(den))
outname = 'full_res_%s_%s_t%04d_%s'%(simname, subset, frame, resolution_name)
fig.savefig(outname)
