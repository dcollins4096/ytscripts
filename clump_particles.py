import pdb
#from yt.mods import *
from yt.visualization.plot_modifications import *
from yt.visualization import _MPL
import pyximport; pyximport.install()
import particle_ops
import particle_clump_mask
import numpy as na
from yt.visualization.plot_modifications import PlotCallback
#reload( particle_clump_mask )

y_dict = [2,2,1]
x_dict = [1,0,0]
debug = 0
print "CLOWN reload 3"
def particles_from_clump(clump):
    out_ind = na.array([],dtype='int64')
    out_x = na.array([])
    out_y = na.array([])
    out_z = na.array([])
    for grid, mask in clump.blocks:

        xpos = grid['particle_position_x']
        ypos = grid['particle_position_y']
        zpos = grid['particle_position_z']
        #pdb.set_trace()

        grid_left = grid.LeftEdge
        grid_dx = grid.dds
        this_cut_mask = mask.astype('int32')

        particle_mask = particle_clump_mask.particle_clump_mask_go( 
            xpos , ypos , zpos, grid_left, grid_dx 
            , this_cut_mask
        )
        particle_mask = particle_mask == 1
        out_ind = na.append(out_ind,grid['particle_index'][particle_mask].astype('int64'))
        out_x =   na.append(out_x,grid['particle_position_x'][particle_mask] )
        out_y =   na.append(out_y,grid['particle_position_y'][particle_mask] )
        out_z =   na.append(out_z,grid['particle_position_z'][particle_mask])

    return out_ind, out_x, out_y, out_z



class DaveParticleCallback(PlotCallback):
    _type_name = "dave_particles"
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




