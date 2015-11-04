# These are the van Genuchten (1980) equations
# The input is matric potential, psi and the hydraulic parameters.
# psi must be sent in as a numpy array.
# The pars variable is like a MATLAB structure.


def thetaFun(psi,pars):
  Se=(1+abs(psi*pars.alpha)**pars.n)**(-pars.m)
  Se[psi>=0]=1.
  return pars.thetaR+(pars.thetaS-pars.thetaR)*Se
  
def CFun(psi,pars):
  Se=(1+abs(psi*pars.alpha)**pars.n)**(-pars.m)
  Se[psi>=0]=1.
  dSedh=pars.alpha*pars.m/(1-pars.m)*Se**(1/pars.m)*(1-Se**(1/pars.m))**pars.m
  # dSedh(psi>=0)=0;
  return Se*pars.Ss+(pars.thetaS-pars.thetaR)*dSedh
  
def KFun(psi,pars):
  Se=(1+abs(psi*pars.alpha)**pars.n)**(-pars.m)
  Se[psi>=0]=1.
  return pars.Ks*Se**pars.neta*(1-(1-Se**(1/pars.m))**pars.m)**2
  
def setpars():
  class Args: pass
  pars=Args()
  pars.thetaR=float(raw_input("thetaR = "))
  pars.thetaS=float(raw_input("thetaS = "))
  pars.alpha=float(raw_input("alpha = "))
  pars.n=float(raw_input("n = "))
  pars.m=1-1/pars.n
  pars.Ks=float(raw_input("Ks = "))
  pars.neta=float(raw_input("neta = "))
  pars.Ss=float(raw_input("Ss = "))
  return pars
  
def PlotProps(pars):
  import numpy as np
  import pylab as pl
  import vanGenuchten as vg
  psi=np.linspace(-10,2,200)
  pl.figure
  pl.subplot(3,1,1)
  pl.plot(psi,vg.thetaFun(psi,pars))
  pl.ylabel(r'$\theta(\psi) [-]$')
  pl.subplot(3,1,2)
  pl.plot(psi,vg.CFun(psi,pars))
  pl.ylabel(r'$C(\psi) [1/m]$')
  pl.subplot(3,1,3)
  pl.plot(psi,vg.KFun(psi,pars))
  pl.xlabel(r'$\psi [m]$')
  pl.ylabel(r'$K(\psi) [m/d]$')
  #pl.show()
  
def HygieneSandstone():
  class Args: pass
  pars=Args()
  pars.thetaR=0.153
  pars.thetaS=0.25
  pars.alpha=0.79
  pars.n=10.4
  pars.m=1-1/pars.n
  pars.Ks=1.08
  pars.neta=0.5
  pars.Ss=0.000001
  return pars
  
def TouchetSiltLoam():
  class Args: pass
  pars=Args()
  pars.thetaR=0.19
  pars.thetaS=0.469
  pars.alpha=0.5
  pars.n=7.09
  pars.m=1-1/pars.n
  pars.Ks=3.03
  pars.neta=0.5
  pars.Ss=0.000001
  return pars

def SiltLoamGE3():
  class Args: pass
  pars=Args()
  pars.thetaR=0.131
  pars.thetaS=0.396
  pars.alpha=0.423
  pars.n=2.06
  pars.m=1-1/pars.n
  pars.Ks=0.0496
  pars.neta=0.5
  pars.Ss=0.000001
  return pars
  
def GuelphLoamDrying():
  class Args: pass
  pars=Args()
  pars.thetaR=0.218
  pars.thetaS=0.520
  pars.alpha=1.15
  pars.n=2.03
  pars.m=1-1/pars.n
  pars.Ks=0.316
  pars.neta=0.5
  pars.Ss=0.000001
  return pars
  
def GuelphLoamWetting():
  class Args: pass
  pars=Args()
  pars.thetaR=0.218
  pars.thetaS=0.434
  pars.alpha=2.0
  pars.n=2.76
  pars.m=1-1/pars.n
  pars.Ks=0.316
  pars.neta=0.5
  pars.Ss=0.000001
  return pars
  
def BeitNetofaClay():
  class Args: pass
  pars=Args()
  pars.thetaR=0.
  pars.thetaS=0.446
  pars.alpha=0.152
  pars.n=1.17
  pars.m=1-1/pars.n
  pars.Ks=0.00082
  pars.neta=0.5
  pars.Ss=0.000001
  return pars