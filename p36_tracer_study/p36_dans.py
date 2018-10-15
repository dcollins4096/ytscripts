from matplotlib import use; use('Agg')
#from matplotlib.backends.backend_agg import FigureCanvasAgg 
import yt
import matplotlib.pyplot as plt
import numpy as np
#---------------------------------------- define tracer fields------------------#
def tracer_number_density(field,data):
    norm = data.get_field_parameter('particle_norm')
    cic = data[('deposit','all_cic')].in_units('code_density')
    number_density = cic/norm
    return number_density
these_units='1' #/code_length**3'
yt.add_field('tracer_number_density',function=tracer_number_density, #units=these_units,
validators=[yt.ValidateParameter("particle_norm")])
#-------------------------------------------------------------------------------#


#---------------------------------------- define tracer fields------------------#
def gas_number_density(field,data):
    norm = data.get_field_parameter('particle_norm')
    cic = data['density'].in_units('code_density')
    g_number_density = cic/norm
    return g_number_density
these_units='1/code_length**3'
yt.add_field('gas_number_density',function=gas_number_density, units=these_units,
validators=[yt.ValidateParameter("particle_norm")])
#-------------------------------------------------------------------------------#

if 0:
    dans_dir = '/scratch1/dkl14b/histo_plot_2'
    ds = yt.load('%s/DD0073/data0073'%(dans_dir))
    ds1 = yt.load('%s/DD0073_MC/data0073'%(dans_dir))
    ds2 = yt.load('%s/DD0073_Flux/data0073'%(dans_dir))
frame = 73
ds  = yt.load('/scratch1/dcollins/Paper19/u07-64-driving-nograv/DD%04d/data%04d'%(frame,frame))
ds1 = yt.load('/scratch1/dcollins/Paper19/u12_montecarlo/DD%04d/data%04d'%(frame,frame))
ds2 = yt.load('/scratch1/dcollins/Paper19/u11_u07_newtracer/DD%04d/data%04d'%(frame,frame))

ad = ds.all_data()
ad1 = ds1.all_data()
ad2 = ds2.all_data()
my_norm = (8*ds.index.grids[0]['particle_mass'][0]/ds.index.grids[0].dds.prod()).in_units('code_density')
ad.set_field_parameter('particle_norm',my_norm)
ad1.set_field_parameter('particle_norm',my_norm)
ad2.set_field_parameter('particle_norm',my_norm)


# for default method
prof_tracer = yt.create_profile(ad,'tracer_number_density',fields='cell_volume',weight_field=None)
#prof_tracer.set_x_unit('1/(code_length**3)')
plot_specs = [dict(color='r')]


# for Monte Carlo method
prof_tracer1 = yt.create_profile(ad1,'tracer_number_density',fields='cell_volume',weight_field=None)
#prof_tracer1.set_x_unit('1/(code_length**3)')
plot_specs = [dict(color='b')]

# for the Flux method
prof_tracer2 = yt.create_profile(ad2,'tracer_number_density',fields='cell_volume',weight_field=None)
#prof_tracer2.set_x_unit('1/(code_length**3)')
plot_specs = [dict(color='k')]


#prof_gas = yt.create_profile(ad,'gas_number_density',fields='cell_volume',weight_field=None)
prof_gas = yt.create_profile(ad,'density',fields='cell_volume',weight_field=None)
#prof_gas.set_x_unit('1/(code_length**3)')
plot_specs.append(dict(color='g'))
       # plot = yt.ProfilePlot.from_profiles([prof_tracer,prof_gas],plot_specs=plot_specs)
plt.clf()
if 0:
    b = max(prof_tracer['cell_volume'])
    a = 0.5*(prof_tracer.x_bins[1:]+prof_tracer.x_bins[0:-1])
    for i in range(0,63):
        if prof_tracer['cell_volume'][i] == b:
            scaler = 1/a[i]

    b1 = max(prof_tracer1['cell_volume'])
    a1 = 0.5*(prof_tracer1.x_bins[1:]+prof_tracer1.x_bins[0:-1])
    for i in range(0,63):
            if prof_tracer1['cell_volume'][i] == b1:
                    scaler1 = 1/a1[i]

    b2 = max(prof_tracer2['cell_volume'])
    a2 = 0.5*(prof_tracer2.x_bins[1:]+prof_tracer2.x_bins[0:-1])
    for i in range(0,63):
            if prof_tracer2['cell_volume'][i] == b2:
                    scaler2 = 1/a2[i]

    bg = max(prof_gas['cell_volume'])
    ag = 0.5*(prof_gas.x_bins[1:]+prof_gas.x_bins[0:-1])
    for i in range(0,63):
            if prof_gas['cell_volume'][i] == bg:
                    scalerg = 1/ag[i]
print ds.current_time
print ds1.current_time
print ds2.current_time


scaler = 1.0
scaler1=1.0
scaler2=1.0
scalerg=1.0
 # ploting default
plt.plot(0.5*(prof_tracer.x_bins[1:]+prof_tracer.x_bins[0:-1])*scaler,prof_tracer['cell_volume'],c='r',label = None)

# ploting MC
plt.plot(0.5*(prof_tracer1.x_bins[1:]+prof_tracer1.x_bins[0:-1])*scaler1,prof_tracer1['cell_volume'],c='b',label = None)

# ploting flux
plt.plot(0.5*(prof_tracer2.x_bins[1:]+prof_tracer2.x_bins[0:-1])*scaler2,prof_tracer2['cell_volume'],c='k',label = None)


# ploting gas
plt.plot(0.5*(prof_gas.x_bins[1:]+prof_gas.x_bins[0:-1])*scalerg,prof_gas['cell_volume'],c='g',
                 label = None)
plt.xscale('log')
plt.yscale('log')
#plt.legend(loc=1)
plt.legend(['default','MC','Flux ave.','gas'],loc=0)
plt.ylabel(r'$V(\rho)$')
plt.xlabel(r'CIC Number density')
fname   = "histo_t_%04d_Num_density_daves.pdf"%(frame)
print fname
print plt.savefig(fname)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------#
if 0:
    ds3 = yt.load('DD0050/data0050')
    ds4 = yt.load('DD0050_MC/data0050')
    ds5 = yt.load('DD0050_Flux/data0050')

    ad3 = ds3.all_data()
    ad4 = ds4.all_data()
    ad5 = ds5.all_data()
    ad3.set_field_parameter('particle_norm',
       8*ds3.index.grids[0]['particle_mass'][0])
    ad4.set_field_parameter('particle_norm',
       8*ds4.index.grids[0]['particle_mass'][0])
    ad5.set_field_parameter('particle_norm',
       8*ds5.index.grids[0]['particle_mass'][0])



# for default method
    prof_tracer3 = yt.create_profile(ad3,'tracer_number_density',fields='cell_volume',weight_field=None)
    prof_tracer3.set_x_unit('1/(code_length**3)')
    plot_specs = [dict(color='r')]


# for Monte Carlo method
    prof_tracer4 = yt.create_profile(ad4,'tracer_number_density',fields='cell_volume',weight_field=None)
    prof_tracer4.set_x_unit('1/(code_length**3)')
    plot_specs = [dict(color='b')]

# for the Flux method
    prof_tracer5 = yt.create_profile(ad5,'tracer_number_density',fields='cell_volume',weight_field=None)
    prof_tracer5.set_x_unit('1/(code_length**3)')
    plot_specs = [dict(color='k')]


    prof_gas3 = yt.create_profile(ad3,'gas_number_density',fields='cell_volume',weight_field=None)
    prof_gas3.set_x_unit('1/(code_length**3)')
    plot_specs.append(dict(color='g'))
           # plot = yt.ProfilePlot.from_profiles([prof_tracer,prof_gas],plot_specs=plot_specs)
    plt.clf()
    b3 = max(prof_tracer3['cell_volume'])
    a3 = 0.5*(prof_tracer3.x_bins[1:]+prof_tracer3.x_bins[0:-1])
    for i in range(0,63):
            if prof_tracer3['cell_volume'][i] == b3:
                    scaler3 = 1/a3[i]

    b4 = max(prof_tracer4['cell_volume'])
    a4 = 0.5*(prof_tracer4.x_bins[1:]+prof_tracer4.x_bins[0:-1])
    for i in range(0,63):
            if prof_tracer4['cell_volume'][i] == b4:
                    scaler4 = 1/a4[i]

    b5 = max(prof_tracer5['cell_volume'])
    a5 = 0.5*(prof_tracer5.x_bins[1:]+prof_tracer5.x_bins[0:-1])
    for i in range(0,63):
            if prof_tracer5['cell_volume'][i] == b5:
                    scaler5 = 1/a5[i]

    bg3 = max(prof_gas3['cell_volume'])
    ag3 = 0.5*(prof_gas3.x_bins[1:]+prof_gas3.x_bins[0:-1])
    for i in range(0,63):
            if prof_gas3['cell_volume'][i] == bg3:
                    scalerg3 = 1/ag3[i]
# ploting default
    plt.plot(0.5*(prof_tracer3.x_bins[1:]+prof_tracer3.x_bins[0:-1])*scaler3,prof_tracer3['cell_volume'],c='r',label = None)

# ploting MC
    plt.plot(0.5*(prof_tracer4.x_bins[1:]+prof_tracer4.x_bins[0:-1])*scaler4,prof_tracer4['cell_volume'],c='b',label = None)

# ploting flux
    plt.plot(0.5*(prof_tracer5.x_bins[1:]+prof_tracer5.x_bins[0:-1])*scaler5,prof_tracer5['cell_volume'],c='k',label = None)

# ploting gas
    plt.plot(0.5*(prof_gas3.x_bins[1:]+prof_gas3.x_bins[0:-1])*scalerg3,prof_gas3['cell_volume'],c='g',
                     label = None)

    print ds3.current_time
    print ds4.current_time
    print ds5.current_time

    plt.xscale('log')
    plt.yscale('log')
#plt.legend(loc=1)
    plt.ylabel(r'$V(\rho)$')
    plt.xlabel(r'CIC number density')
    plt.legend(['default','MC','Flux ave.','gas'])
    fname   = "histo_t_0050_Num_density.pdf"
    print fname
    print plt.savefig(fname)

