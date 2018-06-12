import pdb
import numpy as np 
nar=np.array
import enzo_write
reload(enzo_write)
class waves():
    def rotate(self,k_all, wave=None):
        if wave is None:
            wave = self.wave
        rank = 3
        B0 = self['b0']
        B0_mag = (B0**2).sum()**0.5
        B_unit = B0/B0_mag
        knorm = (k_all[0,...]**2+k_all[1,...]**2+k_all[2,...]**2)**0.5
        ok=knorm>0

        #compute unit vectors
        a_unit = np.zeros_like(k_all)
        
        for dim in range(rank):
            a_unit[dim,...][ok] = k_all[dim,...][ok]/knorm[ok]

        b_unit = np.zeros_like(k_all)
        B0_dot_a = B_unit[0]*a_unit[0,...]+B_unit[1]*a_unit[1,...]+B_unit[2]*a_unit[2,...]
        for dim in range(rank):
            b_unit[dim,...] = -B0_dot_a*a_unit[dim,...]+B_unit[dim]


        b_norm = (b_unit[0,...]**2+b_unit[1,...]**2+b_unit[2,...]**2)**0.5
        ok=b_norm>0
        b_unit_tmp = np.zeros_like(b_unit)
        for dim in range(rank):
            b_unit[dim,...][ok] = b_unit[dim,...][ok]/ b_norm[ok]


        c_unit = np.zeros_like(k_all)
        c_unit[0,...] = a_unit[1,...]*b_unit[2,...]-a_unit[2,...]*b_unit[1,...]
        c_unit[1,...] = a_unit[2,...]*b_unit[0,...]-a_unit[0,...]*b_unit[2,...]
        c_unit[2,...] = a_unit[0,...]*b_unit[1,...]-a_unit[1,...]*b_unit[0,...]

        #CORRECT FOR k||B0.  Not done yet.
        #along = np.logical_and( np.abs(a_unit[0,...]-b_unit[0]) < dx[0],\
        #                        np.abs(a_unit[1,...]-b_unit[1]) < dx[1])
        #along = np.logical_and( np.abs(a_unit[2,...]-b_unit[2]) < dx[2],\
        #                        along)

        B0hat = np.zeros_like(a_unit)
        B0hat[0,...] = B_unit[0]*a_unit[0,...]+B_unit[1]*a_unit[1,...]+B_unit[2]*a_unit[2,...]
        B0hat[1,...] = B_unit[0]*b_unit[0,...]+B_unit[1]*b_unit[1,...]+B_unit[2]*b_unit[2,...]
        B0hat[2,...] = B_unit[0]*c_unit[0,...]+B_unit[1]*c_unit[1,...]+B_unit[2]*c_unit[2,...]
        self.a_unit=a_unit
        self.b_unit=b_unit
        self.c_unit=c_unit
        self.B0hat=B0hat

#this_system = p49_eigen.waves(bx=1.0,by=0.0,bz=0.5,P=0.6,this_wave='a+')
        scalars={}
        for field in ['d','P','e']:
            scalars[field]=self.quan[field]

        hat_system = waves(bx=B0hat[0,...], by=B0hat[1,...], bz=B0hat[2,...],
                           Gamma=self.Gamma,**scalars)
        
        self.rot={}
        for f in ['d','vx','vy','vz','bx','by','bz','e']:
            self.rot[f]=np.zeros_like(a_unit)

        self.rot['d'] = hat_system.right[wave]['d']
        self.rot['e'] = hat_system.right[wave]['e']
        self.rot['vx']  = a_unit[0,...]*hat_system.right[wave]['vx']
        self.rot['vx'] += b_unit[0,...]*hat_system.right[wave]['vy']
        self.rot['vx'] += c_unit[0,...]*hat_system.right[wave]['vz']

        self.rot['vy']  = a_unit[1,...]*hat_system.right[wave]['vx']
        self.rot['vy'] += b_unit[1,...]*hat_system.right[wave]['vy']
        self.rot['vy'] += c_unit[1,...]*hat_system.right[wave]['vz']

        self.rot['vz']  = a_unit[2,...]*hat_system.right[wave]['vx']
        self.rot['vz'] += b_unit[2,...]*hat_system.right[wave]['vy']
        self.rot['vz'] += c_unit[2,...]*hat_system.right[wave]['vz']

        self.rot['bx']  = a_unit[0,...]*hat_system.right[wave]['bx']
        self.rot['bx'] += b_unit[0,...]*hat_system.right[wave]['by']
        self.rot['bx'] += c_unit[0,...]*hat_system.right[wave]['bz']

        self.rot['by']  = a_unit[1,...]*hat_system.right[wave]['bx']
        self.rot['by'] += b_unit[1,...]*hat_system.right[wave]['by']
        self.rot['by'] += c_unit[1,...]*hat_system.right[wave]['bz']

        self.rot['bz']  = a_unit[2,...]*hat_system.right[wave]['bx']
        self.rot['bz'] += b_unit[2,...]*hat_system.right[wave]['by']
        self.rot['bz'] += c_unit[2,...]*hat_system.right[wave]['bz']
        self.hat_system=hat_system


#print("a_unit",a_unit.flatten())
#print("b_unit",b_unit.flatten())
#print("c_unit",c_unit.flatten())


        """
        plt.clf()
        plt.quiver(k_all[0,:,:,0],k_all[1,:,:,0],a_unit[0,:,:,0],a_unit[1,:,:,0])
        plt.quiver(k_all[0,:,:,0],k_all[1,:,:,0],b_unit[0,:,:,0],b_unit[1,:,:,0])
        plt.savefig('p49_units.png')
        """
#a_unit = khat
#sinTheta b_unit = Bvec - cosTheta a_unit
#b_unit = b_uni/||b_unit||
#c_unit = a_unit cross b_unit
    def __init__(self,this_wave='f+', d=1.0,vx=0.0,vy=0.0,vz=0.0,bx=1.0,by=1.41421,bz=0.5,Gamma=1.6666666667,e=None,P=None, write=True):
        #P = 0.6
        #Gamma=1.6666666667
        self.Gamma=Gamma
        self.wave=this_wave
        self.cubes=None
        if P is not None:
            egas=P/((Gamma-1)*d)
            e = 0.5*(vx*vx+vy*vy+vz*vz)+0.5*(bx*bx+by*by+bz*bz)/d + egas
        self.quan={'d':d,'vx':vx,'vy':vy,'vz':vz,'bx':bx,'by':by,'bz':bz,'P':P,'Gamma':Gamma,'egas':egas,'e':e}
        for a,b in [['d','Density'],['vx','x-velocity'],['vy','y-velocity'],['e','TotalEnergy']]:
            self.quan[b]=self.quan[a]
        self.v0 = np.array([vx,vy,vz])
        self.b0 = np.array([bx,by,bz])
        self.right = eigen(d, vx, vy, vz, bx,by,bz,e, Gamma=Gamma, obj=self)
        self.speeds = speeds(d, vx, vy, vz, bx,by,bz,p=P, Gamma=Gamma)
    def __getitem__(self,field):
        map_to_eigen={'d':0,'vx':2,'vy':3,'vz':4,'bx':-1,'by':5,'bz':6,'e':1}
        this_pert = self.right[self.wave]
        #if hasattr(this_pert['by'],'__getitem__'):
        if field is 'dv':
            return nar([ this_pert['vx'], this_pert['vy'], this_pert['vz']])
        if field is 'db':
            return nar([ 0, this_pert['by'], this_pert['bz']])
        if field is 'de':
            return nar(this_pert['d'])
        if field is 'drho':
            return nar(this_pert['e'])
        if field is 'dall':
            return this_pert
        if field is 'b0':
            return self.b0
        if field is 'v0':
            return self.v0
        if field is 'all':
            return self.quan

    def rot_write(self,k_rot=None,base_size=None,pert=1e-6,directory=".", write=True, wave=None):

        self.rotate(k_rot, wave)
        kint = k_rot.astype(int)


        all_hats = {}
        all_p = {}
        cube_slice=[slice(0,base_size[0]),slice(0,base_size[1]),slice(0,base_size[2])]
        field_list = ['d','vx','vy','vz','bx','by','bz','e']
        for f in field_list:
            all_hats[f]=np.zeros(base_size)*1j
            all_hats[f][kint[0,...],kint[1,...],kint[2,...]] = self.rot[f]*pert
            all_hats[f][-kint[0,...],-kint[1,...],-kint[2,...]] = (all_hats[f][kint[0,...],kint[1,...],kint[2,...]] ).conj()
        for f in  field_list:
            tmp=np.fft.ifftn(all_hats[f])
            real_mean = np.mean(np.abs(tmp.real))
            imag_mean = np.mean(np.abs(tmp.imag))
            if (real_mean+imag_mean)>1e-16:
                if imag_mean/(real_mean+imag_mean) > 1e-9:
                    print("Warning: large imaginary component")
            all_p[f]=tmp.real

        if self.cubes is None:
            self.cubes={}
        map_to_label ={'d':'density','vx':'x-velocity','vy':'y-velocity','vz':'z-velocity',
                'bx':'Bx','by':'By','bz':'Bz','e':'TotalEnergy'}
        face_offset = {'bx':nar([1,0,0]),'by':nar([0,1,0]),'bz':nar([0,0,1])}
        for field in ['d','e','vx','vy','vz','bx','by','bz']:
            size = base_size+face_offset.get(field,0)
            this_set = np.ones(size)*self.quan[field]
            if field in ['vx','vy','vz','e']:
                this_set *= self.cubes['density']
            self.cubes[map_to_label[field]]=this_set
        for field in ['d','e','vx','vy','vz','bx','by','bz']:
            size = base_size+face_offset.get(field,0)
            #if field is not 'bx':
                #the_wave = pert*amplitude*self.right[ map_to_eigen[field] ][wave] 
                #the_wave = pert*amplitude*self.right[self.wave][field]
            self.cubes[map_to_label[field]][cube_slice] += all_p[field]
        for field in ['d','e','vx','vy','vz','bx','by','bz']:
            this_filename = "%s/%s_16.h5"%(directory,map_to_label[field])
            if field in ['vx','vy','vz','e']:
                self.cubes[map_to_label[field]] /= self.cubes[map_to_label['d']]
            if write:
                enzo_write.dump_h5(self.cubes[map_to_label[field]],this_filename)


    def perturb(self,base_size=None,pert=1e-6,wave='a+',directory=".", write=True):

        self.cubes={}
        #map_to_eigen={'d':0,'vx':2,'vy':3,'vz':4,'bx':-1,'by':5,'bz':6,'e':1}
        map_to_label ={'d':'density','vx':'x-velocity','vy':'y-velocity','vz':'z-velocity',
                'bx':'Bx','by':'By','bz':'Bz','e':'TotalEnergy'}
        face_offset = {'bx':nar([1,0,0]),'by':nar([0,1,0]),'bz':nar([0,0,1])}
        for field in ['d','e','vx','vy','vz','bx','by','bz']:
            size = base_size+face_offset.get(field,0)
            this_set = np.ones(size)*self.quan[field]
            if field in ['vx','vy','vz','e']:
                this_set *= self.cubes['density']
            self.cubes[map_to_label[field]]=this_set
        for field in ['d','e','vx','vy','vz','bx','by','bz']:
            size = base_size+face_offset.get(field,0)
            amplitude = np.zeros(size)
            amplitude[:size[0]/2,:,:] = 1
            amplitude[size[0]/2:,:,:] = -1
            if field is not 'bx':
                #the_wave = pert*amplitude*self.right[ map_to_eigen[field] ][wave] 
                the_wave = pert*amplitude*self.right[self.wave][field]
                self.cubes[map_to_label[field]] += the_wave
        for field in ['d','e','vx','vy','vz','bx','by','bz']:
            this_filename = "%s/%s_16.h5"%(directory,map_to_label[field])
            if field in ['vx','vy','vz','e']:
                self.cubes[map_to_label[field]] /= self.cubes[map_to_label['d']]
            if write:
                enzo_write.dump_h5(self.cubes[map_to_label[field]],this_filename)
def eigen(d,vx,vy,vz,bx,by,bz,eng, EquationOfState = 0, Gamma=5./3, IsothermalSoundSpeed=1, obj=None):

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

    betay=np.zeros_like(by)+sqrt2i
    betaz=np.zeros_like(by)+sqrt2i
    alph_f=np.zeros_like(by)+1
    alph_s=np.zeros_like(by)+1
    #non-degenerate points
    non_deg = np.logical_or( by != 0.0 , bz != 0.0 )
    bt =np.zeros_like(by)
    if hasattr(by,'__getitem__'):
        bt[non_deg] = np.sqrt( by[non_deg]*by[non_deg] + bz[non_deg]*bz[non_deg] );
        betay[non_deg] = by[non_deg]/bt[non_deg];
        betaz[non_deg] = bz[non_deg]/bt[non_deg];
        alph_f[non_deg] = np.sqrt( (cf[non_deg]*cf[non_deg]-ca[non_deg]*ca[non_deg])/(cf[non_deg]*cf[non_deg]-cs[non_deg]*cs[non_deg]) );
        alph_s[non_deg] = np.sqrt( (cf[non_deg]*cf[non_deg]-aa[non_deg]*aa[non_deg])/(cf[non_deg]*cf[non_deg]-cs[non_deg]*cs[non_deg]) );
    else:
        if non_deg:
            bt = np.sqrt( by*by + bz*bz );
            betay = by/bt;
            betaz = bz/bt;
            alph_f = np.sqrt( (cf*cf-ca*ca)/(cf*cf-cs*cs) );
            alph_s = np.sqrt( (cf*cf-aa*aa)/(cf*cf-cs*cs) );
  
#    if by == 0.0 and bz == 0.0:
#        betay = sqrt2i
#        betaz = sqrt2i
#        alph_f = 1
#        alph_s = 1
#    else:
#        bt = np.sqrt( by*by + bz*bz );
#        betay = by/bt;
#        betaz = bz/bt;
#        alph_f = np.sqrt( (cf*cf-ca*ca)/(cf*cf-cs*cs) );
#        alph_s = np.sqrt( (cf*cf-aa*aa)/(cf*cf-cs*cs) );
  
    #the vectors
    right={}
    for w in ['f-','a-','s-','c','f+','a+','s+']:
        right[w] = {}
        for f in ['d','vx','vy','vz','bx','by','bz','e']:
            right[w][f]=np.zeros_like(by)
    map_to_eigen={'d':0,'vx':2,'vy':3,'vz':4,'bx':-1,'by':5,'bz':6,'e':1}
  
    #fast, left
    right['f-']['d'] = alph_f
    if EquationOfState == 0 :
        right['f-']['e'] =  alph_f*0.5*(vx*vx+vy*vy+vz*vz) 
        right['f-']['e'] += alph_f*cf*cf*og1 - alph_f*cf*vx + alph_s*ca*sbx*(betay*vy + betaz*vz) 
        right['f-']['e'] += (Gamma-2)*og1*alph_f*(cf*cf-aa*aa)
    right['f-']['vx'] = alph_f*(vx - cf);
    right['f-']['vy'] = alph_f*vy + alph_s*betay*ca*sbx;
    right['f-']['vz'] = alph_f*vz + alph_s*betaz*ca*sbx;
    right['f-']['by'] = alph_s*betay*cf*sqrtDi;
    right['f-']['bz'] = alph_s*betaz*cf*sqrtDi;


    #alfven][left
    #right['a-']['d'] = 0; already zero
    if EquationOfState == 0:
        right['a-']['e'] = 1*(betaz*vy - betay*vz)*sbx;
    #right['a-']['vx'] = 0; already zero.
    right['a-']['vy'] =  1*betaz*sbx;
    right['a-']['vz'] = -1*betay*sbx;
    right['a-']['by'] = -betaz*sqrtDi;
    right['a-']['bz'] = betay*sqrtDi;

    #alfven,right
    #right['a+']['d'] = 0.; already zero
    if EquationOfState == 0.:
        right['a+']['e'] = -1*(betaz*vy - betay*vz)*sbx;
    #right['a+']['vx'] = 0.; already zero.
    right['a+']['vy'] = -1*betaz*sbx;
    right['a+']['vz'] = +1*betay*sbx;
    right['a+']['by'] = betaz*sqrtDi;
    right['a+']['bz'] = -betay*sqrtDi;

  
    #slow,left
    right['s-']['d'] = alph_s;
    if EquationOfState == 0:
        right['s-']['e'] = alph_s*0.5*(vx*vx+vy*vy+vz*vz) + \
            alph_s*cs*cs*og1 - alph_s*cs*vx - alph_f*aa*sbx*(betay*vy + betaz*vz) +\
            (Gamma-2)*og1*alph_s*(cs*cs - aa*aa );
    right['s-']['vx'] = alph_s*(vx-cs);
    right['s-']['vy'] = alph_s*vy - alph_f*betay*aa*sbx;
    right['s-']['vz'] = alph_s*vz - alph_f*betaz*aa*sbx;
    right['s-']['by'] = -alph_f*betay*aa*aa*sqrtDi/cf;
    right['s-']['bz'] = -alph_f*betaz*aa*aa*sqrtDi/cf;
  
    #entropy (no entropy wave in isothermal MHD.)(Or hydro,for that matter)
    if EquationOfState == 0:
        right['c']['d'] = 1;
        right['c']['e'] = 0.5*(vx*vx+vy*vy+vz*vz);
        right['c']['vx'] = vx;
        right['c']['vy'] = vy;
        right['c']['vz'] = vz;
        right['c']['by'] = 0;
        right['c']['bz'] = 0;
  
    #slow,right
    right['s+']['d'] = alph_s;
    if EquationOfState  == 0:
        right['s+']['e'] = alph_s*0.5*(vx*vx+vy*vy+vz*vz) + alph_s*cs*cs*og1 + \
            alph_s*cs*vx + alph_f*aa*sbx*(betay*vy + betaz*vz) + (Gamma-2)*og1*alph_s*(cs*cs - aa*aa );
    right['s+']['vx'] = alph_s*(vx+cs);
    right['s+']['vy'] = alph_s*vy + alph_f*betay*aa*sbx;
    right['s+']['vz'] = alph_s*vz + alph_f*betaz*aa*sbx;
    right['s+']['by'] = -alph_f*betay*aa*aa*sqrtDi/cf;
    right['s+']['bz'] = -alph_f*betaz*aa*aa*sqrtDi/cf;
  


    #fast, right
    right['f+']['d'] = alph_f;
    if EquationOfState  == 0:
        right['f+']['e'] = alph_f*0.5*(vx*vx+vy*vy+vz*vz) + \
        alph_f*cf*cf*og1 + alph_f*cf*vx - alph_s*ca*sbx*(betay*vy + betaz*vz) + \
        (Gamma-2)*og1*alph_f*(cf*cf-aa*aa);
    right['f+']['vx'] = alph_f*(vx + cf);
    right['f+']['vy'] = alph_f*vy - alph_s*betay*ca*sbx;
    right['f+']['vz'] = alph_f*vz - alph_s*betaz*ca*sbx;
    right['f+']['by'] = alph_s*betay*cf*sqrtDi;
    right['f+']['bz'] = alph_s*betaz*cf*sqrtDi;
    return right
def speeds(d,vx,vy,vz,bx,by,bz,eng=None,p=None, EquationOfState = 0, Gamma=5./3, IsothermalSoundSpeed=1):

    #float d, float vx, float vy, float vz, 
    #float bx, float by, float bz, float eng, 
    #float right[][7])

    #Eigen vectors taken from Ryu & Jones, ApJ 442:228-258, 1995
    #Normalization starting on p. 231.

    bp = 0.5*(bx*bx + by*by + bz*bz );
    sqrt2 = np.sqrt(2.0);
    sqrt2i= 1.0/sqrt2;
    sqrtD = np.sqrt(d);
    sqrtDi = 1.0/sqrtD;
    sbx  = np.sign(bx);
    og1  = 1.0/(Gamma - 1);
    if p is not None:
        Egas=p/((Gamma-1)*d)
        eng = 0.5*(vx*vx+vy*vy+vz*vz)+0.5*(bx*bx+by*by+bz*bz)/d + Egas


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


