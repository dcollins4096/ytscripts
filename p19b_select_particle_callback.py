import yt
from yt.visualization.plot_modifications import *
import pyximport; pyximport.install()
import particle_ops
import numpy as na
dbg = 0
def deposit_target_particles_1(field, data):
    #pdb.set_trace()
    indices_late = data.get_field_parameter('indices_late')
    blank = np.zeros(data.ActiveDimensions, dtype='float64')
    if not hasattr(indices_late,'sum'):
        return blank
    if data.NumberOfParticles == 0: return blank
    if dbg > 0:
        if hasattr(indices_late,'sum'):
            print "indices inside", indices_late.sum()
        else:
            print "indices inside", indices_late
            return blank
    #int 32 because it's just a flag.
    #mask_to_get = np.zeros(indices_late.shape, dtype='int32')
    mask_to_get = data.get_field_parameter('mask_to_get')
    t0 = data.get_field_parameter('timer')
    my_indices = data['particle_index'].astype('int64')
    #print "  min thing"
    #mask_to_get[ indices_late < my_indices.min()] = 1
    #mask_to_get[ indices_late > my_indices.max()] = 1
    #print "  left", mask_to_get.size - mask_to_get.sum()
    #print "  search"
    found_any, mask = particle_ops.mask_particles(
        indices_late, my_indices, mask_to_get)
    data.set_field_parameter('mask_to_get',mask_to_get)
    pos_to_get = data['particle_position'][mask == 1]
    d = data.deposit(pos_to_get, method = "count")
    d = data.ds.arr(d, input_units = "cm**-3")
    t1 = time.time() 
    #print "  DT", t1-t0[-1]
    data.set_field_parameter('timer',t0+[t1])
    return data.apply_units(d, field.units)

yt.add_field(      ("deposit","deposit_target_particles"),
#yt.add_field(      'dp1',
         function = deposit_target_particles_1,
         validators = [yt.ValidateSpatial(), 
                       yt.ValidateParameter('indices_late'), 
                       yt.ValidateParameter('mask_to_get'), 
                       yt.ValidateGridType()],
         display_name = "target_particles")


class SelectParticleCallback(PlotCallback):
    _type_name = "select_particles"
    region = None
    _descriptor = None
    def __init__(self, width, p_size=1.0, col='k', marker='o', stride=1.0,
                 ptype='all', stars_only=False, dm_only=False,
                 minimum_mass=None,bool=None, indices=None, xyz=None):
        """
        Adds particle positions, based on a thick slab along *axis* with a
        *width* along the line of sight.  *p_size* controls the number of
        pixels per particle, and *col* governs the color.  *ptype* will
        restrict plotted particles to only those that are of a given type.
        *minimum_mass* will require that the particles be of a given mass,
        calculated via ParticleMassMsun, to be plotted.
        """
        PlotCallback.__init__(self)
        self.width = width
        self.p_size = p_size
        self.color = col
        self.marker = marker
        self.stride = stride
        self.ptype = ptype
        self.stars_only = stars_only
        self.dm_only = dm_only
        self.minimum_mass = minimum_mass
        self.bool = bool
        self.indices = indices
        self.xyz = xyz

    def __call__(self, plot):
        data = plot.data
        if iterable(self.width):
            self.width = np.float64(plot.data.ds.quan(self.width[0], self.width[1]))
        # we construct a recantangular prism
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        if self.bool is None:
            reg = self._get_region((x0,x1), (y0,y1), plot.data.axis, data)
        else:
            reg = self.bool
        ax = data.axis

        xax = plot.data.ds.coordinates.x_axis[ax]
        yax = plot.data.ds.coordinates.y_axis[ax]
        axis_names = plot.data.ds.coordinates.axis_name
        field_x = "particle_position_%s" % axis_names[xax]
        field_y = "particle_position_%s" % axis_names[yax]
        pt = self.ptype
        gg = ( ( reg[pt, field_x] >= x0 ) & ( reg[pt, field_x] <= x1 )
           &   ( reg[pt, field_y] >= y0 ) & ( reg[pt, field_y] <= y1 ) )
        if self.indices != None:
            print "particle callback warning: this might get expensive"
            mask_to_get = na.zeros(self.indices.shape, dtype='int32')
            found_any, mask = particle_ops.mask_particles(
                self.indices.astype('int64'), reg['particle_index'].astype('int64'), mask_to_get)
            gg = ( gg & (mask == 1) )

        print "nparticles do it.", mask_to_get.sum()
        if False:
            if self.ptype is not None:
                gg &= (reg["particle_type"] == self.ptype)
                if gg.sum() == 0: return
            if self.stars_only:
                gg &= (reg["creation_time"] > 0.0)
                if gg.sum() == 0: return
            if self.dm_only:
                gg &= (reg["creation_time"] <= 0.0)
                if gg.sum() == 0: return
            if self.minimum_mass is not None:
                gg &= (reg["ParticleMassMsun"] >= self.minimum_mass)
                if gg.sum() == 0: return
        plot._axes.hold(True)
        px, py = self.convert_to_plot(plot,
                    [np.array(reg[pt, field_x][gg][::self.stride]),
                     np.array(reg[pt, field_y][gg][::self.stride])])
        plot._axes.scatter(px, py, edgecolors='None', marker=self.marker,
                           s=self.p_size, c=self.color)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

    def _get_region(self, xlim, ylim, axis, data):
        LE, RE = [None]*3, [None]*3
        ds=data.ds
        xax = ds.coordinates.x_axis[axis]
        yax = ds.coordinates.y_axis[axis]
        zax = axis
        LE[xax], RE[xax] = xlim
        LE[yax], RE[yax] = ylim
        LE[zax] = data.center[zax].ndarray_view()- self.width*0.5
        RE[zax] = data.center[zax].ndarray_view()+ self.width*0.5
        if self.region is not None \
            and na.all(self.region.left_edge <= LE) \
            and na.all(self.region.right_edge >= RE):
            return self.region
        self.region = data.ds.region( data.center, LE, RE)
        return self.region




