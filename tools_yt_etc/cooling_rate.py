########################################################################
#
# Cooling rate example script
#
#
# Copyright (c) 2013-2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from matplotlib import pyplot
import numpy as np
import os
import yt

from pygrackle import \
    chemistry_data, \
    setup_fluid_container

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc

def make_chem(grackle_data_file,with_radiative_cooling=0,primordial_chemistry=3,metal_cooling=1,UVbackground=1,
              self_shielding_method=0,H2_self_shielding=0,comoving_coordinates=0,a_units=1.0,a_value=1.0,current_redshift = 0.):
    print("making chem.")

    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = with_radiative_cooling
    my_chemistry.primordial_chemistry = primordial_chemistry
    my_chemistry.metal_cooling = metal_cooling
    my_chemistry.UVbackground = UVbackground
    my_chemistry.self_shielding_method = self_shielding_method
    my_chemistry.H2_self_shielding =H2_self_shielding 
    #grackle_dir = os.path.dirname(os.path.dirname(os.path.dirname(
    #    os.path.dirname(os.path.abspath(__file__)))))
    #my_chemistry.grackle_data_file = os.sep.join(
    #    [grackle_dir, "input", "CloudyData_UVB=HM2012.h5"])
    my_chemistry.grackle_data_file = grackle_data_file

    # Set units
    #current_redshift = 0.
    my_chemistry.comoving_coordinates = comoving_coordinates # proper units
    my_chemistry.a_units = a_units
    my_chemistry.a_value = a_value / (1.0 + current_redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
    my_chemistry.time_units = sec_per_Myr          # 1 Gyr in s
    my_chemistry.velocity_units = my_chemistry.a_units * \
        (my_chemistry.length_units / my_chemistry.a_value) / \
        my_chemistry.time_units
    print(my_chemistry.grackle_data_file)
    return my_chemistry

def plot_stuff(my_chemistry,density=1.6737352238051868e-24, plot_object=pyplot, **plot_args):
    print(my_chemistry.grackle_data_file)
    # Call convenience function for setting up a fluid container.
    # This container holds the solver parameters, units, and fields.
    temperature = np.logspace(1, 9, 200)
    fc = setup_fluid_container(my_chemistry,
                               temperature=temperature,density=density,
                               converge=True)

    print(my_chemistry.grackle_data_file)
    fc.calculate_temperature()
    fc.calculate_cooling_time()

    print(my_chemistry.grackle_data_file)
    density_proper = fc["density"] / \
        (my_chemistry.a_units *
         my_chemistry.a_value)**(3*my_chemistry.comoving_coordinates)
    cooling_rate = fc.cooling_units * fc["energy"] / \
        np.abs(fc["cooling_time"]) / density_proper

    data = {}
    t_sort = np.argsort(fc["temperature"])
    for field in fc.density_fields:
        data[field] = yt.YTArray(fc[field][t_sort] *
                                 my_chemistry.density_units, "g/cm**3")
    data["energy"]       = yt.YTArray(fc["energy"][t_sort], "erg/g")
    data["temperature"]  = yt.YTArray(fc["temperature"][t_sort], "K")
    data["pressure"]     = yt.YTArray(fc["pressure"][t_sort], "dyne/cm**2")
    data["cooling_time"] = yt.YTArray(fc["cooling_time"][t_sort], "s")
    data["cooling_rate"] = yt.YTArray(cooling_rate[t_sort], "erg*cm**3/s")

    plot_object.loglog(data["temperature"], data["cooling_rate"],**plot_args)
    plot_object.xlabel('T [K]')
    plot_object.ylabel('$\\Lambda$ [erg s$^{-1}$ cm$^{3}$]')

    # save data arrays as a yt dataset
    if 'PRIMORDIAL_CHEM' in os.environ:
        ds_name = 'cooling_rate.pc%s.h5' % os.environ['PRIMORDIAL_CHEM']
        im_name = 'cooling_rate.pc%s.png' % os.environ['PRIMORDIAL_CHEM']
    else:
        ds_name = 'cooling_rate.h5'
        im_name = 'cooling_rate.png'
    plot_object.savefig(im_name)
    print(im_name)
    return data
    #yt.save_as_dataset({}, ds_name, data)
