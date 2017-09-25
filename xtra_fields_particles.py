

from yt.data_objects.particle_filters import add_particle_filter
def formed_star(pfilter, data):
    filter = data["all", "creation_time"] > 0
    return filter

add_particle_filter("formed_star", function=formed_star, filtered_type='all',
                                        requires=["creation_time"])
