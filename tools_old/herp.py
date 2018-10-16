
    ds = fct * (data['dx'].flat[0])**2.0
    f  =         data["MagneticPressure"][ip1 ,ip1 ,ip1 ]
    f +=         data["MagneticPressure"][ip1 ,ip1 ,i   ]
    f +=         data["MagneticPressure"][ip1 ,ip1 ,im1 ]

    f +=         data["MagneticPressure"][ip1 ,i   ,ip1 ]
    f +=         data["MagneticPressure"][ip1 ,i   ,i   ]
    f +=         data["MagneticPressure"][ip1 ,i   ,im1 ]

    f +=         data["MagneticPressure"][ip1 ,im1 ,ip1 ]
    f +=         data["MagneticPressure"][ip1 ,im1 ,i   ]
    f +=         data["MagneticPressure"][ip1 ,im1 ,im1 ]

    f +=         data["MagneticPressure"][i   ,ip1 ,ip1 ]
    f +=         data["MagneticPressure"][i   ,ip1 ,i   ]
    f +=         data["MagneticPressure"][i   ,ip1 ,im1 ]

    f +=         data["MagneticPressure"][i   ,i   ,ip1 ]
    #f +=         data["MagneticPressure"][i   ,i   ,i   ]  all at once 
    f +=         data["MagneticPressure"][i   ,i   ,im1 ]

    f +=         data["MagneticPressure"][i   ,im1 ,ip1 ]
    f +=         data["MagneticPressure"][i   ,im1 ,i   ]
    f +=         data["MagneticPressure"][i   ,im1 ,im1 ]

    f +=         data["MagneticPressure"][im1 ,ip1 ,ip1 ]
    f +=         data["MagneticPressure"][im1 ,ip1 ,i   ]
    f +=         data["MagneticPressure"][im1 ,ip1 ,im1 ]

    f +=         data["MagneticPressure"][im1 ,i   ,ip1 ]
    f +=         data["MagneticPressure"][im1 ,i   ,i   ]
    f +=         data["MagneticPressure"][im1 ,i   ,im1 ]

    f +=         data["MagneticPressure"][im1 ,im1 ,ip1 ]
    f +=         data["MagneticPressure"][im1 ,im1 ,i   ]
    f +=         data["MagneticPressure"][im1 ,im1 ,im1 ]

    f -=      26*data["MagneticPressure"][i   ,i   ,i   ] 
