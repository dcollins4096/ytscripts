if 1:
        a =( (ds.domain_right_edge-ds.domain_left_edge)/ds.domain_dimensions ).max()
        L = ds.domain_left_edge+1.5*a
        R = ds.domain_right_edge - 1.5*a
        C = 0.5*(L+R)
        reg = ds.region( C,L,R)

