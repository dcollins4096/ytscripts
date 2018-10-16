"""
Get energies for a structure.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.utilities.data_point_utilities import FindBindingEnergy

    
def energies_codeunits(clump, G, use_thermal_energy=True, use_particles=False, truncate=False):
    """puts kinetic, thermal, gravitational, and magnetic energies in code units into a dictionary"""

    energy_dict = {}
    if False:
        use_particles &= \
          ("all", "particle_mass") in clump.data.ds.field_info
    
    bulk_velocity = clump.quantities.bulk_velocity(use_particles=use_particles)
    energy_dict['bulk_velocity'] = bulk_velocity

    kinetic = 0.5 * (clump["gas", "cell_mass"] *
        ((bulk_velocity[0] - clump["gas", "velocity_x"])**2 +
         (bulk_velocity[1] - clump["gas", "velocity_y"])**2 +
         (bulk_velocity[2] - clump["gas", "velocity_z"])**2)).sum()

    energy_dict['kinetic'] = kinetic
    if use_thermal_energy:
        if clump.ds['EquationOfState'] == 0:
            thermal =  (clump["gas", "cell_mass"] * clump["gas", "thermal_energy"]).sum()
        else:
            thermal = 1.5*(clump["cell_volume"]*clump["density"]).sum()
    energy_dict['thermal'] = thermal

    if use_particles:
        kinetic += 0.5 * (clump["all", "particle_mass"] *
            ((bulk_velocity[0] - clump["all", "particle_velocity_x"])**2 +
             (bulk_velocity[1] - clump["all", "particle_velocity_y"])**2 +
             (bulk_velocity[2] - clump["all", "particle_velocity_z"])**2)).sum()

    potential = G * FindBindingEnergy(clump["gas", "cell_mass"],
                          clump["index", "x"],
                          clump["index", "y"],
                          clump["index", "z"],
                          truncate, (kinetic / G))


    
    if use_particles:
        potential += clump.data.ds.quan(G *
            FindBindingEnergy(
                clump["all", "particle_mass"],
                clump["all", "particle_position_x"],
                clump["all", "particle_position_y"],
                clump["all", "particle_position_z"],
                truncate, ((kinetic - potential) / G)), kinetic)

    energy_dict['gravitational'] = potential


    energy_dict['bulk_velocity'] = bulk_velocity
    return energy_dict
