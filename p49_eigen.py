import numpy as np 
def speeds(d,vx,vy,vz,bx,by,bz,eng, EquationOfState = 0, Gamma=5./3, IsothermalSoundSpeed=1):

    #float d, float vx, float vy, float vz, 
    #float bx, float by, float bz, float eng, 
    #float right[][7])

    #Eigen vectors taken from Ryu & Jones, ApJ 442:228-258, 1995
    #Normalization starting on p. 231.

    p = -1
    bp = 0.5*(bx*bx + by*by + bz*bz );
    sqrt2 = np.sqrt(2.0);
    sqrt2i= 1.0/sqrt2;
    sqrtD = np.sqrt(d);
    sqrtDi = 1.0/sqrtD;
    sbx  = np.sign(bx);
    og1  = 1.0/(Gamma - 1);

    if EquationOfState == 0 :
        p = (Gamma -1.0 ) * (eng - 0.5* d * (vx*vx + vy*vy + vz * vz ) - 0.5*(bx*bx + by*by + bz*bz));
    else:
        p = IsothermalSoundSpeed*IsothermalSoundSpeed * d;

    #compute wave speeds
    if EquationOfState == 0:
        aa = np.sqrt( 1.0*Gamma* p/d )
    else:
        aa = IsothermalSoundSpeed;

    cs =np.sqrt( 0.5*( aa*aa + 2*bp/d -np.sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );
    ca =np.sqrt( bx*bx/d ); 
    cf =np.sqrt( 0.5*( aa*aa + 2*bp/d +np.sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );
    return {'cf':cf,'cs':cs,'ca':ca,'aa':aa}
def eigen(d,vx,vy,vz,bx,by,bz,eng, EquationOfState = 0, Gamma=5./3, IsothermalSoundSpeed=1):

    #float d, float vx, float vy, float vz, 
    #float bx, float by, float bz, float eng, 
    #float right[][7])

    #Eigen vectors taken from Ryu & Jones, ApJ 442:228-258, 1995
    #Normalization starting on p. 231.

    p = -1
    bp = 0.5*(bx*bx + by*by + bz*bz );
    sqrt2 = np.sqrt(2.0);
    sqrt2i= 1.0/sqrt2;
    sqrtD = np.sqrt(d);
    sqrtDi = 1.0/sqrtD;
    sbx  = np.sign(bx);
    og1  = 1.0/(Gamma - 1);

    if EquationOfState == 0 :
        p = (Gamma -1.0 ) * (eng - 0.5* d * (vx*vx + vy*vy + vz * vz ) - 0.5*(bx*bx + by*by + bz*bz));
    else:
        p = IsothermalSoundSpeed*IsothermalSoundSpeed * d;

    #compute wave speeds
    if EquationOfState == 0:
        aa = np.sqrt( 1.0*Gamma* p/d )
    else:
        aa = IsothermalSoundSpeed;

    cs =np.sqrt( 0.5*( aa*aa + 2*bp/d -np.sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );
    ca =np.sqrt( bx*bx/d ); 
    cf =np.sqrt( 0.5*( aa*aa + 2*bp/d +np.sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );
    #compute ancilary values
    #The normalization of alph_f may change. This normalization uses Ryu
    #& Jones, but Balsara may be more robust.

    if by == 0.0 and bz == 0.0:
        betay = sqrt2i
        betaz = sqrt2i
        alph_f = 1
        alph_s = 1
    else:
        bt = np.sqrt( by*by + bz*bz );
        betay = by/bt;
        betaz = bz/bt;
        alph_f = np.sqrt( (cf*cf-ca*ca)/(cf*cf-cs*cs) );
        alph_s = np.sqrt( (cf*cf-aa*aa)/(cf*cf-cs*cs) );
  
    #the vectors
  
    #fast, left
    right=np.zeros([7,7])
    right[0,0] = alph_f
    if EquationOfState == 0 :
        right[1,0] =  alph_f*0.5*(vx*vx+vy*vy+vz*vz) 
        right[1,0] += alph_f*cf*cf*og1 - alph_f*cf*vx + alph_s*ca*sbx*(betay*vy + betaz*vz) 
        right[1,0] += (Gamma-2)*og1*alph_f*(cf*cf-aa*aa)
    right[2,0] = alph_f*(vx - cf);
    right[3,0] = alph_f*vy + alph_s*betay*ca*sbx;
    right[4,0] = alph_f*vz + alph_s*betaz*ca*sbx;
    right[5,0] = alph_s*betay*cf*sqrtDi;
    right[6,0] = alph_s*betaz*cf*sqrtDi;

    #alfven][left
    right[0,1] = 0;
    if EquationOfState == 0:
        right[1,1] = 1*(betaz*vy - betay*vz)*sbx;
    right[2,1] = 0;
    right[3,1] =  1*betaz*sbx;
    right[4,1] = -1*betay*sbx;
    right[5,1] = betaz*sqrtDi;
    right[6,1] = -betay*sqrtDi;

  
    #slow,left
    right[0,2] = alph_s;
    if EquationOfState == 0:
        right[1,2] = alph_s*0.5*(vx*vx+vy*vy+vz*vz) + \
            alph_s*cs*cs*og1 - alph_s*cs*vx - alph_f*aa*sbx*(betay*vy + betaz*vz) +\
            (Gamma-2)*og1*alph_s*(cs*cs - aa*aa );
    right[2,2] = alph_s*(vx-cs);
    right[3,2] = alph_s*vy - alph_f*betay*aa*sbx;
    right[4,2] = alph_s*vz - alph_f*betaz*aa*sbx;
    right[5,2] = -alph_f*betay*aa*aa*sqrtDi/cf;
    right[6,2] = -alph_f*betaz*aa*aa*sqrtDi/cf;
  
    #entropy (no entropy wave in isothermal MHD.)(Or hydro,for that matter)
    if EquationOfState == 0:
        right[0,3] = 1;
        right[1,3] = 0.5*(vx*vx+vy*vy+vz*vz);
        right[2,3] = vx;
        right[3,3] = vy;
        right[4,3] = vz;
        right[5,3] = 0;
        right[6,3] = 0;
  
    #slow,right
    right[0,4] = alph_s;
    if EquationOfState  == 0:
        right[1,4] = alph_s*0.5*(vx*vx+vy*vy+vz*vz) + alph_s*cs*cs*og1 + \
            alph_s*cs*vx + alph_f*aa*sbx*(betay*vy + betaz*vz) + (Gamma-2)*og1*alph_s*(cs*cs - aa*aa );
    right[2,4] = alph_s*(vx+cs);
    right[3,4] = alph_s*vy + alph_f*betay*aa*sbx;
    right[4,4] = alph_s*vz + alph_f*betaz*aa*sbx;
    right[5,4] = -alph_f*betay*aa*aa*sqrtDi/cf;
    right[6,4] = -alph_f*betaz*aa*aa*sqrtDi/cf;
  
    #alfven,right
    right[0,5] = 0;
    if EquationOfState == 0:
        right[1,5] = -1*(betaz*vy - betay*vz)*sbx;
    right[2,5] = 0;
    right[3,5] = -1*betaz*sbx;
    right[4,5] = +1*betay*sbx;
    right[5,5] = betaz*sqrtDi;
    right[6,5] = -betay*sqrtDi;


    #fast, right
    right[0,6] = alph_f;
    if EquationOfState  == 0:
        right[1,6] = alph_f*0.5*(vx*vx+vy*vy+vz*vz) + \
        alph_f*cf*cf*og1 + alph_f*cf*vx - alph_s*ca*sbx*(betay*vy + betaz*vz) + \
        (Gamma-2)*og1*alph_f*(cf*cf-aa*aa);
    right[2,6] = alph_f*(vx + cf);
    right[3,6] = alph_f*vy - alph_s*betay*ca*sbx;
    right[4,6] = alph_f*vz - alph_s*betaz*ca*sbx;
    right[5,6] = alph_s*betay*cf*sqrtDi;
    right[6,6] = alph_s*betaz*cf*sqrtDi;
    return right
