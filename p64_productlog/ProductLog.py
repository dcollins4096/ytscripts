from go import *
from scipy.special import lambertw  
plt.clf()
x = np.arange(-1,1,0.01)
"""
1/Sqrt[2 \[Pi] \[Sigma]^2] (\[Alpha]/\[Epsilon] ProductLog[
     0, \[Epsilon]/(c \[Rho]1)]/(
    1 + Log[\[Epsilon]/(c \[Rho]1)/
      ProductLog[0, \[Epsilon]/(c \[Rho]1)]])
     Exp[-((\[Alpha] Log[\[Epsilon]/(c \[Rho]0)/
          ProductLog[0, \[Epsilon]/(c \[Rho]1)]] - \[Mu])^2/(
      2 \[Sigma]^2))] - \[Alpha]/\[Epsilon] ProductLog[-1, \
\[Epsilon]/(c \[Rho]1)]/(
    1 + Log[\[Epsilon]/(c \[Rho]1)/
      ProductLog[-1, \[Epsilon]/(c \[Rho]1)]])
     Exp[-((\[Alpha] Log[\[Epsilon]/(c \[Rho]0)/
          ProductLog[-1, \[Epsilon]/(c \[Rho]1)]] - \[Mu])^2/(
      2 \[Sigma]^2))] HeavisideTheta[-\[Epsilon]])
      """

def ProductLog(k,x):
    return lambertw(x,k)
def HeavisideTheta(x):
    return np.heaviside(x,0)

def g0(epsilon=1, mu=1, sigma=1,alpha=1,c=1,rho0=1,rho1=1):
    out = 1/np.sqrt(2* np.pi*sigma**2)*(alpha/epsilon* ProductLog( 0, epsilon/(c*rho1))/( 1 + np.log(epsilon/(c*rho1)/ ProductLog(0, epsilon/(c*rho1))))*\
                                       np.exp(-((alpha*np.log(epsilon/(c*rho0)/ ProductLog(0, epsilon/(c*rho1))) - mu)**2/( 2*sigma**2)))
                                       - alpha/epsilon*ProductLog(-1, epsilon/(c*rho1))/( 1 + np.log(epsilon/(c*rho1)/ ProductLog(-1, epsilon/(c*rho1))))*\
                                       np.exp(-((alpha*np.log(epsilon/(c*rho0)/ ProductLog(-1, epsilon/(c*rho1))) - mu)**2/( 2*sigma**2)))* HeavisideTheta(-epsilon))
    return out
def g1(epsilon=1, mu=1, sigma=1,alpha=1,c=1,rho0=1,rho1=1):
    e1 = epsilon/(c*rho1)
    e0 = epsilon/(c*rho0)
    out = 1/np.sqrt(2* np.pi*sigma**2)*(alpha/epsilon* ProductLog( 0, e1)/( 1 + np.log(e1/ ProductLog(0, e1)))*\
                                       np.exp(-((alpha*np.log(e0/ ProductLog(0, e1)) - mu)**2/( 2*sigma**2)))
                                       - alpha/epsilon*ProductLog(-1, e1)/( 1 + np.log(e1/ ProductLog(-1, e1)))*\
                                       np.exp(-((alpha*np.log(e0/ ProductLog(-1, e1)) - mu)**2/( 2*sigma**2)))* HeavisideTheta(-epsilon))
    return out

def g2(epsilon=1, mu=1, sigma=1,alpha=1,c=1,rho0=1,rho1=1):
    e1 = epsilon/(c*rho1)
    e0 = epsilon/(c*rho0)
    out = 1/np.sqrt(2* np.pi*sigma**2)*(alpha/epsilon* ProductLog( 0, e1)/( 1 + np.log(e1/ ProductLog(0, e1)))*\
                                       np.exp(-((alpha*np.log(e0/ ProductLog(0, e1)) - mu)**2/( 2*sigma**2)))
                                       - alpha/epsilon*ProductLog(-1, e1)/( 1 + np.log(e1/ ProductLog(-1, e1)))*\
                                       np.exp(-((alpha*np.log(e0/ ProductLog(-1, e1)) - mu)**2/( 2*sigma**2)))* HeavisideTheta(-epsilon))
    return out
plt.clf()
x=np.arange(-1,6,0.1)
plt.plot( x, g1( x,mu=1, sigma=0.25,c=1,rho0=1,rho1=1))
plt.plot( x, g1( x,mu=1, sigma=0.5,c=1,rho0=1,rho1=1))
plt.plot( x, g1( x,mu=1, sigma=1.0,c=1,rho0=1,rho1=1))
plt.plot( x, g1( x,mu=1, sigma=2.0,c=1,rho0=1,rho1=1))
plt.savefig('%s/RESEARCH3/Paper64_AcousticEnergy/p64_test2.png')
