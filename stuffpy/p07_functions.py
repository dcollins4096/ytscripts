def imshow_symlog(my_matrix, vmin, vmax, logthresh=5):
    img=imshow( my_matrix ,
                vmin=float(vmin), vmax=float(vmax),
                norm=matplotlib.colors.SymLogNorm(10**-logthresh) )

    maxlog=int(np.ceil( np.log10(vmax) ))
    minlog=int(np.ceil( np.log10(-vmin) ))
    #generate logarithmic ticks 
    tick_locations=([-(10**x) for x in xrange(minlog,-logthresh-1,-1)]
                    +[0.0]
                    +[(10**x) for x in xrange(-logthresh,maxlog+1)] )

    cb=colorbar(ticks=tick_locations)
    return img,cb
# FUNCTIONS

def _Q(field,data):

        '''This function calculates the Stokes Parameter "Q" by using magnetic
        field components along the line of sight. For exemple: line of sight is
        along the x-axis, therefore, the mag-field components used are y and z.

        Makes use of the polarization factor "epsilon" using a power exponent.
        Epsilon is defined within the function.
        '''

        p = data.get_field_parameter('p')
        n0 = data.get_field_parameter('n0')
        if p is None:
            p=-1
        if n0 is None:
            n0 = 1

        epsilon = np.ones(data['density'].shape)
        n = data['density']
        epsilon[ n > n0 ] = ((n/n0)**p)[ n > n0 ]

        return ( epsilon*n * ( ((data['magnetic_field_'+d2])**2.0) - \
                             ((data['magnetic_field_'+d3])**2.0) ) )

yt.add_field('Q', units='g**2/(cm**4*s**2)', function=_Q)

################################################################################

def _U(field,data):

        '''This function calculates the Stokes Parameter "U" by using magnetic
        field components along the line of sight. For exemple: line of sight is
        along the x-axis, therefore, the mag-field components used are y and z.

        Makes use of the polarization factor "epsilon" using a power exponent.
        Epsilon is defined within the function.
        '''

        p = data.get_field_parameter('p')
        n0 = data.get_field_parameter('n0')
        if p is None:
            p=-1
        if n0 is None:
            n0 = 1

        epsilon = np.ones(data['density'].shape)
        n = data['density']
        epsilon[ n > n0 ] = ((n/n0)**p)[ n > n0 ]


        # I reversed y,z to be z,y to correspond to the dimensional
        # coord(x,z,y) that appears to be set for the data.
        return ( 2.0 * epsilon * n * (data['magnetic_field_'+d2]) * \
                                 (data['magnetic_field_'+d3]) )

yt.add_field('U', units='g**2/(cm**4*s**2)', function=_U)

#theta = 0.5 arctan( u/q)

def segments(histo,n_segments,cumul=False):
    """Spits out indices where the histogram has *n_segments* equal count segments.
    *cumul* if its already a cumulative histogram."""
    if not cumul:
        cumu_hist = np.cumsum(histo)
        cumu_hist = cumu_hist/float(cumu_hist.max())
    else:
        cumu_hist = histo/histo.max()
    frac_boundaries = np.arange(1,n_segments)/float(n_segments)
    indices= [np.where( (cumu_hist-frac) >= 0.0)[0][0] for frac in frac_boundaries]
    #pdb.set_trace()
    return indices, cumu_hist

