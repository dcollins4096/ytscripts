from go import *

def addit(exp,signed_k, label,v,phase_angle=0,verbose=True):
    try:
        k = tuple(np.abs(signed_k))
        phase_sign = np.sign(signed_k[-1])
    except:
        k = np.abs(signed_k)
        phase_sign = np.sign(signed_k)
    if phase_sign==0: phase_sign=1
    if verbose:
        exp['labs'] = exp.get('labs',{})
        exp['full_history']=exp.get('full_history',{})

    if (np.abs(v) > 1e-10).any():
        this_value = v*np.exp(1j*phase_angle*phase_sign)
        if (np.abs(k)<1e-16).all():
            mult=2 
            this_value *= mult
            this_value=this_value.real
        exp[k] = exp.get(k,0)+this_value
        if verbose:
            exp['labs'][label]=k
            exp['full_history'][label] = v*np.exp(1j*phase_angle*phase_sign)

def modes_1(KA=0,KB=0,KC=0,A0=0,A1=0,B0=0,B1=0,C0=0,C1=0, Atheta=0,Btheta=0,Ctheta=0):


    expect={}

    addit( expect, np.zeros_like(KA),       "0",                0.5* A0*B0*C0)
    addit( expect, KA,      "KA",        C0*B0*(0.5*A1), Atheta)
    addit( expect, KB,      "KB",        A0*C0*(0.5*B1), Btheta)
                                      
    addit( expect, KA+KB,   "KA+KB",     C0*(0.5*A1)*(0.5*B1), Atheta+Btheta)
    addit( expect, KA-KB,   "KA-KB",     C0*(0.5*A1)*(0.5*B1), Atheta-Btheta)
    addit( expect, KC,      "KC",        A0*B0*(0.5*C1), Ctheta)
    addit( expect, KB+KC,   "KB+KC",     A0*(0.5*B1)*(0.5*C1), Btheta+Ctheta)
    addit( expect, KB-KC,   "KB-KC",     A0*(0.5*B1)*(0.5*C1), Btheta-Ctheta)
                                      
    addit( expect, KA-KC,   "KA-KC",     B0*(0.5*A1)*(0.5*C1), (Atheta-Ctheta))
    addit( expect, KA+KC,   "KA+KC",     B0*(0.5*A1)*(0.5*C1), Atheta+Ctheta)
                                      
    addit( expect, KA+KB+KC,"KA+KB+KC",  (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta+Btheta+Ctheta)
    addit( expect, KA-KB+KC,"KA-KB+KC",  (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta-Btheta+Ctheta)
    addit( expect, -KA+KB+KC,"KA-KB-KC", (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta-Btheta-Ctheta)
    addit( expect, KA+KB-KC,"KA+KB-KC",  (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta+Btheta-Ctheta)
    return expect

def modes_2(A=[],B=[],K=[],Theta=[],verbose=True):
    # signal = prod(i)(Ai + Bi exp(1j ki x) )
    nmodes = len(A)
    imode = np.arange(nmodes)
    #A0 = A[0]
    #B0 = A[1]
    #A1 = B[0]
    #B1 = B[1]
    #KA = K[0]
    #KB = K[1]
    #C0 = A[2]
    #C1 = B[2]
    #KC = K[2]
    #Atheta= Theta[0]
    #Btheta= Theta[1]
    #Ctheta= Theta[2]
    A=nar(A)
    B=nar(B)
    K=nar(K)
    Theta=nar(Theta)
    def mode_label(*args_in): 
        alpha = 'ABCDE'
        args = list(args_in)
        output = "K%s"%alpha[args.pop(0)]
        pm_d={1:'+',-1:'-'}
        while len(args):
            output += pm_d[ args.pop(0)]
            output += "K%s"%alpha[args.pop(0)]

        return output
    def these_parts(A,*args):
        ok = np.ones_like(A,dtype='bool')
        for n in args:
            ok = np.logical_and( ok, imode != n)
        return np.prod( A[ok])

    expect={}

    addit( expect, np.zeros_like(K[0]),       "0",                0.5* np.prod(A),verbose=verbose)

    #addit( expect, KA,      "KA",        C0*B0*(0.5*A1), Atheta)
    for n in range(nmodes):
        mode = K[n]
        theta= Theta[n]
        amp = np.prod( A[ imode != n])*B[n]*0.5
        addit( expect, mode,      mode_label(n),        amp,theta,verbose=verbose)
                                      

    #def these_parts(A,n,m):
    #    ok = np.logical_and( imode != n, imode != m)
    #    return np.prod( A[ok])

    for n in range(nmodes):
        for m in range(n+1,nmodes):
            for pm in [1,-1]: 
                mode = K[n] + pm*K[m]
                theta = Theta[n] + pm*Theta[m]
                amp = these_parts(A,n,m)*(0.5*B[n])*(0.5*B[m])
                addit( expect, mode,   mode_label(n,pm,m),    amp, theta,verbose=verbose)

    #addit( expect, KA+KB,   "KA+KB",     C0*(0.5*A1)*(0.5*B1), Atheta+Btheta)
    #addit( expect, KA-KB,   "KA-KB",     C0*(0.5*A1)*(0.5*B1), Atheta-Btheta)
    #addit( expect, KB+KC,   "KB+KC",     A0*(0.5*B1)*(0.5*C1), Btheta+Ctheta)
    #addit( expect, KB-KC,   "KB-KC",     A0*(0.5*B1)*(0.5*C1), Btheta-Ctheta)
    #addit( expect, KA-KC,   "KA-KC",     B0*(0.5*A1)*(0.5*C1), (Atheta-Ctheta))
    #addit( expect, KA+KC,   "KA+KC",     B0*(0.5*A1)*(0.5*C1), Atheta+Ctheta)


    if nmodes >= 3:
        for n in range(nmodes):
            for m in range(n+1,nmodes):
                for ell in range(m+1,nmodes):
                    for pm1, pm2 in [ [1,1], [-1,1],[1,-1],[-1,-1]]:
                        mode = K[n] + pm1*K[m] + pm2*K[ell]
                        theta = Theta[n] + pm1*Theta[m] + pm2*Theta[ell]
                        amp = these_parts(A,n,m,ell)*(0.5*B[n])*(0.5*B[m])*(0.5*B[ell])
                        addit( expect, mode, mode_label(n,pm1,m,pm2,ell) ,    amp, theta,verbose=verbose)


                  
    if nmodes >= 4:
        plus_minus = []
        for pm1 in [1,-1]:
            for pm2 in [1,-1]:
                for pm3 in [1,-1]:
                    plus_minus.append( [pm1,pm2,pm3])
        for n in range(nmodes):
            for m in range(n+1,nmodes):
                for ell in range(m+1,nmodes):
                    for oo in range(ell+1,nmodes):
                        for pm1, pm2, pm3 in plus_minus:
                            mode = K[n] + pm1*K[m] + pm2*K[ell]+pm3*K[oo]
                            theta = Theta[n] + pm1*Theta[m] + pm2*Theta[ell]+pm3*Theta[oo]
                            amp = these_parts(A,n,m,ell,oo)*(0.5*B[n])*(0.5*B[m])*(0.5*B[ell])*(0.5*(B[oo]))
                            addit( expect, mode, mode_label(n,pm1,m,pm2,ell,pm3,oo) ,    amp, theta,verbose=verbose)
    if nmodes >= 5:
        plus_minus = []
        for pm1 in [1,-1]:
            for pm2 in [1,-1]:
                for pm3 in [1,-1]:
                    for pm4 in [1,-1]:
                        plus_minus.append( [pm1,pm2,pm3,pm4])
        for n1 in range(nmodes):
            for n2 in range(n1+1,nmodes):
                for n3 in range(n2+1,nmodes):
                    for n4 in range(n3+1,nmodes):
                        for n5 in range(n4+1,nmodes):
                            for pm1, pm2, pm3,pm4 in plus_minus:
                                mode = K[n1] + pm1*K[n2] + pm2*K[n3]+pm3*K[n4]+pm4*K[n5]
                                theta = Theta[n1] + pm1*Theta[n2] + pm2*Theta[n3]+pm3*Theta[n4]+pm4*Theta[n5]
                                amp = these_parts(A,n1,n2,n3,n4,n5)*(0.5*B[n1])*(0.5*B[n2])*(0.5*B[n3])*(0.5*(B[n4]))*(0.5*(B[n5]))
                                addit( expect, mode, mode_label(n1,pm1,n2,pm2,n3,pm3,n4,pm4,n5) ,    amp, theta,verbose=verbose)
    #addit( expect, KA+KB-KC,"KA+KB-KC",  (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta+Btheta-Ctheta)



                                      
    #addit( expect, KA+KB+KC,"KA+KB+KC",  (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta+Btheta+Ctheta)
    #addit( expect, KA-KB+KC,"KA-KB+KC",  (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta-Btheta+Ctheta)
    #addit( expect, -KA+KB+KC,"-KA+KB+KC", (0.5*A1)*(0.5*B1)*(0.5*C1), -Atheta+Btheta+Ctheta)
    #addit( expect, KA+KB-KC,"KA+KB-KC",  (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta+Btheta-Ctheta)
    return expect
