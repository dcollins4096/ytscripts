if 'ef' not in dir():
    execfile('go')
from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter, RadMC3DSource
#
# Need to make the thing centered.
#

basedir = '/Users/dcollins/scratch/RADMC/'
if 1:
    simname = '%s/StarParticles/plrd01000'%basedir
    sim = 'StarParticles'
    position_cm = [0.0, 0.0, 0.0]
if 0:
    sim = 'StupidSphere'
    frame = 1
    simname = '%s/StupidSphere/DD%04d/data%04d'%(basedir,frame,frame)
    position_cm = [0.5, 0.5, 0.5]
if 0:
    sim = 'IsolatedGalaxy'
    frame = 30
    simname = '%s/IsolatedGalaxy/galaxy%04d/galaxy%04d'%(basedir,frame,frame)
    position_cm = [0.0, 0.0, 0.0]

if 1:
    ds = yt.load(simname)

    dust_to_gas = 0.01
    def _DustDensity(field, data):
        return dust_to_gas * data["density"]
    ds.add_field(("gas", "dust_density"), function=_DustDensity, units="g/cm**3")
    ad=ds.all_data()
    profile = yt.create_profile(ad, 'dust_density', 'cell_mass', weight_field=None)
    plt.clf()
    xb = 0.5*(profile.x_bins[1:]+profile.x_bins[0:-1])
    cell_mass = profile['cell_mass']
    plt.plot(xb,cell_mass)
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel('dust deinsity'); plt.ylabel('cell_mass')
    plt.savefig('dust_pdf'+sim)

    writer = RadMC3DWriter(ds)
    writer.write_amr_grid()
    writer.write_dust_file(("gas", "dust_density"), "dust_density.inp")
    radius_cm = 6.96e10
    mass_g = 1.989e33
    temperature_K = 5780.0
    star = RadMC3DSource(radius_cm, mass_g, position_cm, temperature_K)
    sources_list = [star]
    wavelengths_micron = np.logspace(-1.0, 4.0, 1000)
    writer.write_source_files(sources_list, wavelengths_micron)

else:
    plt.clf()
    from yt.analysis_modules.radmc3d_export.api import read_radmc3d_image
    header, image = read_radmc3d_image("%s/%s/image.out"%(basedir,sim))

    Nx = header['Nx']
    Ny = header['Ny']

    if 0:
        x_hi = 0.5*header["pixel_size_cm_x"]*Nx
        x_lo = -x_hi
        y_hi = 0.5*header["pixel_size_cm_y"]*Ny
        y_lo = -y_hi
    if 1:
        x_hi = header["pixel_size_cm_x"]*Nx
        x_lo = 0.0
        y_hi = header["pixel_size_cm_y"]*Ny
        y_lo = 0.0

    X = np.linspace(x_lo, x_hi, Nx)
    Y = np.linspace(y_lo, y_hi, Ny)

    zeros = image==0
    out = np.log10(image)
    out[zeros] = np.log10(image[image>0].min())
    plt.imshow(out,interpolation='nearest',origin='lower')

    #plt.pcolormesh(X, Y, out ) #, cmap='hot')
    #plt.pcolormesh(X, Y, image, cmap='hot')
    cbar = plt.colorbar()
    plt.axis((x_lo, x_hi, y_lo, y_hi))
    ax = plt.gca()
    ax.set_xlabel(r"$x$ (cm)")
    ax.set_ylabel(r"$y$ (cm)")
    cbar.set_label(r"Log Intensity (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$)")
    plt.savefig('dust_continuum_%s.png'%sim)
