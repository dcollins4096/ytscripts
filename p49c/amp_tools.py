from go import *

def addit(exp,signed_k, label,v,phase_angle=0):
    try:
        k = tuple(np.abs(signed_k))
        phase_sign = np.sign(signed_k[-1])
    except:
        k = np.abs(signed_k)
        phase_sign = np.sign(signed_k)
    if phase_sign==0: phase_sign=1
    exp['labs'] = exp.get('labs',{})
    exp['full_history']=exp.get('full_history',{})

    if np.abs(v) > 1e-10:
        this_value = v*np.exp(1j*phase_angle*phase_sign)
        if np.abs(k)<1e-16:
            mult=2 
            this_value *= mult
            this_value=this_value.real
        exp[k] = exp.get(k,0)+this_value
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
    addit( expect, -KA+KB+KC,"-KA+KB+KC", (0.5*A1)*(0.5*B1)*(0.5*C1), -Atheta+Btheta+Ctheta)
    addit( expect, KA+KB-KC,"KA+KB-KC",  (0.5*A1)*(0.5*B1)*(0.5*C1), Atheta+Btheta-Ctheta)
    return expect

