
execfile('go')
frame = 60; basedir = '/scratch1/dcollins/Paper08/B02/512/'
fname = '%s/RS%04d/restart%04d'%(basedir, frame,frame)
peaks = fPickle.load('p12_b02_512_n0060_peaks2.pickle')

ds = yt.load(fname)
width = 0.05
field='velocity_divergence'
ax=0
peak = peaks[0]
data = ds.sphere(peak,(width,'code_length'))
proj = ds.proj(field, ax, data_source = data, center = peak)
pw = proj.to_pw(center = peak, width = 2*width)
pw.save('p12_oklahoma_fix.png')
