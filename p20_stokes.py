
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
        p = 0
    if n0 is None:
        n0 = 1
    n0 = 1

    epsilon = np.ones(data['density'].shape)
    n = data['density']
    epsilon[ n > n0 ] = ((n.v/n0)**p)[ n > n0 ]
    return ( epsilon*n * ( ((data['Bx'])**2.0) - \
                                 ((data['By'])**2.0) ) )
#    return ( ( ((data['Bx'])**2.0) - \
#                                 ((data['By'])**2.0) ) )


yt.add_field('Q', units='g**2/(cm**4*s**2)', function=_Q)
#yt.add_field('Q', units='code_magnetic**2', function=_Q)

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
        p = 0
    if n0 is None:
        n0 = 1

    n0 = 1

    epsilon = np.ones(data['density'].shape)
    n = data['density']
    epsilon[ n > n0 ] = ((n.v/n0)**p)[ n > n0 ]


# I reversed y,z to be z,y to correspond to the dimensional
# coord(x,z,y) that appears to be set for the data.
   return ( 2.0 * epsilon * n * (data['Bx']) * \
           (data['By']) )
#    return ( 2.0 * (data['Bx']) * \
#            (data['By']) )

yt.add_field('U', units='g**2/(cm**4*s**2)', function=_U)
#yt.add_field('U', units='code_magnetic**2', function=_U)
#       dproj.set_field_parameter('n0',ds.quan(1.0,'g/cm**3'))
#       dproj.set_field_parameter('p',-1.0)

