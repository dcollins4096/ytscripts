
import yt
from yt.visualization.plot_modifications import *
import pyximport; pyximport.install()
import particle_ops
import time
import numpy as na

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
        #plot._axes.hold(True)
        px, py = self.convert_to_plot(plot,
                    [np.array(reg[pt, field_x][gg][::self.stride]),
                     np.array(reg[pt, field_y][gg][::self.stride])])
        plot._axes.scatter(px, py, edgecolors='None', marker=self.marker,
                           s=self.p_size, c=self.color)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        #plot._axes.hold(False)

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





class AnotherStreamlineCallback(PlotCallback):
    """
    annotate_streamlines(field_x, field_y, factor=16,
                         density=1, plot_args=None):

    Add streamlines to any plot, using the *field_x* and *field_y*
    from the associated data, skipping every *factor* datapoints like
    'quiver'. *density* is the index of the amount of the streamlines.
    """
    _type_name = "streamlines_dave"
    def __init__(self, field_x, field_y, factor = 16,
                 density = 1, plot_args=None,other_proj=None):
        PlotCallback.__init__(self)
        def_plot_args = {}
        self.field_x = field_x
        self.field_y = field_y
        self.factor = factor
        self.dens = density
        if plot_args is None: plot_args = def_plot_args
        self.plot_args = plot_args
        self.other_proj = other_proj
        
    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        plot._axes.hold(True)
        nx = plot.image._A.shape[0] / self.factor
        ny = plot.image._A.shape[1] / self.factor
        the_plot_data = plot.data
        if self.other_proj is not None:
            the_plot_data = self.other_proj.data_source
        pixX = _MPL.Pixelize(the_plot_data['px'],
                             the_plot_data['py'],
                             the_plot_data['pdx'],
                             the_plot_data['pdy'],
                             the_plot_data[self.field_x],
                             int(nx), int(ny),
                             (x0, x1, y0, y1),).transpose()
        pixY = _MPL.Pixelize(the_plot_data['px'],
                             the_plot_data['py'],
                             the_plot_data['pdx'],
                             the_plot_data['pdy'],
                             the_plot_data[self.field_y],
                             int(nx), int(ny),
                             (x0, x1, y0, y1),).transpose()
        X,Y = (np.linspace(xx0,xx1,nx,endpoint=True),
               np.linspace(yy0,yy1,ny,endpoint=True))
        streamplot_args = {'x': X, 'y': Y, 'u':pixX, 'v': pixY,
                           'density': self.dens}
        streamplot_args.update(self.plot_args)
        plot._axes.streamplot(**streamplot_args)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class DaveParticleCallback_old_and_broken(PlotCallback):
    _type_name = "dave_particles_old_and_broken"
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

