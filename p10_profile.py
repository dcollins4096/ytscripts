
fiducial_dt = 1e-6
def _Qplus(field,data):
    dt = data.get_field_parameter('dt')
    output=data['QInstantaneous']
    negative = output < 0
    output[negative] = 0
    return output

yt.add_field('Qplus', function = _Qplus)

def rainbow_map(n):
    norm = mpl.colors.Normalize()
    norm.autoscale(na.arange(n))
    #cmap = mpl.cm.jet
    color_map = mpl.cm.ScalarMappable(norm=norm,cmap='jet')
    return  color_map.to_rgba

output_template = '/data/astro12_3/student3/runs/%d/DD%04d/data%04d'

all_xbins = []
all_profiles = []
fields = ['Qplus','cell_volume']
for sim in [193]:
    for frame in range(0,1100,100):
        fname = output_template%(sim,frame,frame)
        ds = yt.load(fname)
        reg = ds.all_data()
        weight_field = None
        prof = yt.create_profile(reg,fields[0],fields[1],weight_field=weight_field)
        all_xbins.append(prof.x_bins)
        all_profiles.append(prof[fields[1]])


plt.clf()
n_prof = len(all_profiles)
rm = rainbow_map(n_prof)
for i in range(n_prof):
    plt.plot( 0.5*(all_xbins[i][1:]+all_xbins[i][0:-1]), all_profiles[i],c=rm(i))

plt.xscale('log'); plt.yscale('log')
plt.xlabel(fields[0]); plt.ylabel(fields[1])
plt.savefig('test.pdf')


