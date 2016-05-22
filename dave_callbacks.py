
from yt.visualization.plot_modifications import *
from yt.visualization import _MPL
from yt.visualization.plot_modifications import PlotCallback

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
