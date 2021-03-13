import numpy as np
from scipy.optimize import root
from scipy.optimize import fsolve
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
from bdcsym import bdcsym
from EOCVWCV import EOCVWCV
from bdcWCV import bdcWCV
from potenfindr import potenfindr
from fzeroforrf import fzeroforrf
from odefuncexfull import odefuncexfull
from odephipf import odephipf
from IntepoAB import IntepoAB
from findBoverC import findBoverC
from othermroot import othermroot
from otherzero import otherzero
import matplotlib.pyplot as plt
from zerointerval import zerointerval
from fzerotx import fzerotx
from rk45 import rk45
from factorfront import factorfront
#from assimulo.solvers import CVode
#from assimulo.solvers import LSODAR
#from assimulo.problem import Explicit_Problem
#import Parameters as P

def PFHYBRID(T_temperature,i_current,L_length,RpH,RpO,PO2EH,PO2EO,x_portion,ENUM,rformer,factortorL,factortorR,formerindex,po2index,IndexInitial,valrtol,valatol,comtol,methodRK):
   
    #global R, F, T, c0, i, iv, idindex, inindex, cL, L
    #PFHYBRID(800,-0.8,0.002,0.1,0.1,1e-23,0.2,0.2,2e+4,0,1,1,1,1,1)
    R=8.314  #J/(mol*K)
    F=96485.0  #C/mol
    Avorgacon=6.02*1e+23 #mol^-1

    #######Operating Conditions of the Cell##########

    T=T_temperature+273.15 #K
    i=i_current #Acm^-2
    eVoverkT=1.0/8.618*1e+5/T
    
    ######## PARAMETERS IN YSZ ######################

    sigmaionYSZ=163*np.exp(-0.79*eVoverkT)*factortorL #OHM^-1cm^-1
    DeYSZ=230.0*np.exp(-2*eVoverkT)*factortorL #cm^2/s
    DhYSZ=0.23*np.exp(-1.15*eVoverkT)*factortorL
    KehYSZ=(1.04*10**23*np.exp(-1.99*eVoverkT))/Avorgacon*1.72*10**21*np.exp(-0.62*eVoverkT)/Avorgacon
    conYSZ=(1.04*10**23*np.exp(-1.99*eVoverkT))/Avorgacon #mol/cm^3
    conYSZh=1.72*10**21*np.exp(-0.62*eVoverkT)/Avorgacon
    KehYSZ=KehYSZ/conYSZ/conYSZ
    AYSZ=5.1385
    CY=2*8/(AYSZ**3*10**(-24))/Avorgacon*0.08/2.16
    #CY=0.0078288; %mol cm-3
    DvYSZ=sigmaionYSZ*R*T/4/F/F/(CY/2)
    #pause
    cGdYSZ=CY #mol/cm3;
    cGdYSZ=cGdYSZ/conYSZ
    bYSZ=cGdYSZ
    #YSZ=P.LeftLayer(T_temperature)
    #print('DeYSZ=',YSZ.De)

    ######## PARAMETERS IN CGO ######################
    muh0=2.4*1e+6*np.exp(-1.45*1.0/8.618*1e+5/1073.15)/1073.15/(1.602*1e-19)/(9.26e+18)/np.exp(-1.08*1.0/8.618*1e+5/1073.15)*1073.15
    muhCGO=muh0/T*np.exp(-1.08*eVoverkT)
    conCGO=4*1e+9*np.exp(-2.51*eVoverkT)/(1.602*1e-19)/(10**(-4301/T+4.1947))/Avorgacon*factortorR #mol/cm^3
    conCGOh=2.4*1e+6*np.exp(-1.45*eVoverkT)/T/muhCGO/(1.602*1e-19)/Avorgacon
    conCGOh=2.4*1e+6*np.exp(-1.45*eVoverkT)/(1392*np.exp(-1.08*eVoverkT))/(1.602*1e-19)/Avorgacon*factortorR
    DvCGO=10**(-3405/T-1.8995)*factortorR
    DeCGO=10**(-4301/T+0.1301)*factortorR #cm^2/s
    DhCGO=muh0*R/F*np.exp(-1.08*eVoverkT)
    #DhCGO=0.12*exp(-1.08*eVoverkT)*factortorR;
    KehCGO=conCGO*conCGOh
    KehCGO=KehCGO/conCGO/conCGO
    #pause
    AGDC=5.418
    cGdCGO=2*8/(AGDC**3*10**(-24))/Avorgacon*0.05/2.0
    #cGdCGO=0.0041789*factortorR; %mol/cm3;
    #conCGO-conCGOh
    cGdCGO=cGdCGO/conCGO
    bCGO=cGdCGO

    ######### parameters on model choosing ##############

    idindex = 1
    inindex = 1

    ######## parameters on the structure of the cell ########
    
    L=L_length #CM
    x0=L*x_portion
    x0=L-x0

    ######## phase field parameter \epsilon #################

    epsilon=L/ENUM

    ####### BOUNDARY CONDITIONS #############################

    if idindex==1 and inindex==1:
       if po2index==0:
          [cce0, cceL]=bdcsym(RpH,RpO,PO2EH,PO2EO,i,R,T,F)
          cce0=np.array([cce0])
          cceL=np.array([cceL])
       elif po2index==1:
            if x_portion>0 and x_portion<1:
               [c0I, c0L]=bdcsym(RpH,RpO,PO2EH,PO2EO,0,R,T,F)
               EOCVH=EOCVWCV(KehYSZ,bYSZ,conYSZ,AYSZ,PO2EH,c0I,R,T,F)
               EOCVO=EOCVWCV(KehCGO,bCGO,conCGO,AGDC,PO2EO,c0L,R,T,F)
               [cce0, cceL]=bdcWCV(RpH,RpO,EOCVH,EOCVO,i,R,T,F)
            elif x_portion>=1:
               [c0I, c0L]=bdcsym(RpH,RpO,PO2EH,PO2EO,0,R,T,F)
               EOCVH=EOCVWCV(KehCGO,bCGO,conCGO,AGDC,PO2EH,c0I,R,T,F)
               EOCVO=EOCVWCV(KehCGO,bCGO,conCGO,AGDC,PO2EO,c0L,R,T,F)
               [cce0, cceL]=bdcWCV(RpH,RpO,EOCVH,EOCVO,i,R,T,F) 
            else:
               [c0I, c0L]=bdcsym(RpH,RpO,PO2EH,PO2EO,0,R,T,F)
               EOCVH=EOCVWCV(KehYSZ,bYSZ,conYSZ,AYSZ,PO2EH,c0I,R,T,F)
               EOCVO=EOCVWCV(KehYSZ,bYSZ,conYSZ,AYSZ,PO2EO,c0L,R,T,F)
               [cce0, cceL]=bdcWCV(RpH,RpO,EOCVH,EOCVO,i,R,T,F)
    else:
       [cce0, cceL]=bdcsym(RpH,RpO,PO2EH,PO2EO,i,R,T,F)
       cce0=np.array([cce0])
       cceL=np.array([cceL])

    c0=cce0
    cL=cceL
    #print('c0=',c0)
    #print('typec0',type(c0))
    #print('shapec0',c0.shape)
    
    Estfactor=factorfront(i_current,L,x_portion,T_temperature)
    po2oeled=np.exp(np.log(PO2EO)-2*R*T/(2*F)*np.arcsinh(i/(2*R*T/2/F/RpO))*4*F/R/T)
    po2guess=po2oeled*np.exp(Estfactor)
    cLguess=np.power(po2guess,(-1/4))
    #print('testinf=',np.sinh(i*RpO*F/R/T))
    #print('Est=',Estfactor)
    #print('po2eled=',po2oeled)
    #print('po2guess=',po2guess)
    #print('CL=',cL)
    #print('clguess=',cLguess)
    ##### Solve r and Ce ################
    #print("test IndexInitial = " + str(IndexInitial))
    if formerindex==1:
       #rinitial=potenfindr(i,RpH,RpO,PO2EH,PO2EO,L*(1-x_portion),T-273.15,1,c0,cLguess)
       rinitial=potenfindr(i,RpH,RpO,PO2EH,PO2EO,L,T-273.15,1,c0,cL)
        #print('rinitial1=', rinitial)
       #rinitial=otherzero(DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, rinitial, 2e+4)
       #print(rinitial)
    elif formerindex==-1:
       rinitial=-0.972
    else:
       rinitial=rformer

    #print(rinitial)
    #print(c0)
    #print(type(c0))
    #rinitial=2.63697
    #[rt, infodict, ier, mesg]=fsolve(fzeroforrf, rinitial, (DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i), fprime=None, full_output=1)
    [lend, rend, riniguess,numberit]=zerointerval(DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol, rinitial,methodRK)
    #print('Lend=', lend)
    #print('rend=', rend)
    if (lend*rend<0) or (numberit>=50):
        xp=0
        po2xp=0
        err=1
        return xp, po2xp, err
    
    #sol=root_scalar(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol), method='bisect', bracket=[lend, rend], xtol=2.2204e-16)
    #r=sol.root
    #r=brentq( lend, rend)
     #   [rt, infodict, ier, mesg]=fsolve(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol), [lend, rend], full_output=1)
    #root=newton(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i), x0=5201.6)
    #if ier==1:
    #   xx=rt.tolist()
    #   r=xx[0]
       #r=sol.root
    
    #print('rbisection=', r)
    
    #sol=root_scalar(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol), method='brentq', bracket=[lend, rend], xtol=2.2204e-16)
    #r=sol.root
    #print('rbrentq=', r)
    #sol=root_scalar(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol), method='brenth', bracket=[lend, rend], xtol=2.2204e-16)
    #r=sol.root
    #print('rbrenth=', r)
    #sol=root_scalar(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol), method='ridder', bracket=[lend, rend], xtol=2.2204e-16)
    #r=sol.root
    #print('rridder=', r)
    
    [rt, infodict, ier, mesg]=fsolve(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol,methodRK), riniguess, full_output=1)
    xx=rt.tolist()
    r=xx[0]
    #sol=root_scalar(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol), method='toms748', bracket=[lend, rend], xtol=2.2204e-16)
    #r=sol.root
    #print('rfsolve=', r)
    r=fzerotx(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol,methodRK), lend, rend,2.2204e-16)
    #print('rfzerotx=',r)
    #r=1.93527
    #sol=root_scalar(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol), method='newton', bracket=None, x0=lend, x1=rend, xtol=2.2204e-16)
    #r=sol.root
    #print('rnewton=', r)
    #sol=root_scalar(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol), method='secant', bracket=None, x0=lend, x1=rend, xtol=2.2204e-16)
    #r=sol.root
    #print('secant=', r)
    #sol=root_scalar(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol), method='halley', bracket=None, x0=lend, x1=rend, xtol=2.2204e-16)
    #r=sol.root
    #print('rhalley=', r)
    #   print('ro=',root)
    #else:
    #    print('Please give another initial guess')
    #    print(mesg)
    #    print(infodict['fvec'])
    #    print(infodict['qtf'])
    #    print(infodict['nfev'])
    #     [lend,rend]=othermroot(DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, rinitial, valrtol, valatol)
    #    print('lend=', lend)
    #     print('rend=', rend)
    #    r=brentq(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol), lend, rend)
        #print('rbrq=', r)
        #[rt, infodict, ier, mesg]=fsolve(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i),  r, full_output=1)
        #xx=rt.tolist()
        #r=xx[0]
        #print('rfsolve=', r)
        
    #[lend,rend]=othermroot(DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, rinitial)
    #r=brentq(lambda tt: fzeroforrf(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i), lend, rend, xtol=1e-12, rtol=1e-12)
    #print('rbrentq=',r)
    #r=5396.01
    #print(lend)
    #print(rend)
        
    if IndexInitial==1:
       y0=cL
    #   xspan=[L:-epsilon/4:0]
       #xspan=np.arange(L, -epsilon/4, -epsilon/4)
       xspan0=np.linspace(0, L, num=4*ENUM+1)
       xspan=L-xspan0
       print('xspan[0]=',xspan[0])
       print('xspan[len]=',xspan[len(xspan)-1])
    elif IndexInitial==-1:
       y0=c0
    #  xspan=[0:epsilon/4:L]
       xspan=np.linspace(0, L, num=4*ENUM+1)

    #print('y0shape=',y0.shape)
    #print('r=',r)

    #####model = Explicit_Problem(lambda x,y: odefuncexfull(x,y,r,DvYSZ, DvCGO, x0, KehYSZ, KehCGO,\
    #####    bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO,epsilon,inindex,\
    #####    idindex, F, i), y0, L)
    #####sim = CVode(model)
    #####sim.rtol=1.e-3
    #####sim.atol=1.e-6
    #####tfinal = 0.0        #Specify the final time
    #####sim.backward = True
    #####xp, Cexp = sim.simulate(tfinal)
    # [xp,Cexp]=ode45(@(x,y) odefuncexfull(r, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon,x,y), xspan, y0);
    #sol=solve_ivp(lambda x,y: odefuncexfull(x,y,r,DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO,epsilon,inindex, idindex, F, i), [xspan[0], xspan[-1]], y0, method='RK45', t_eval=xspan, dense_output=False, events=None, vectorized=False, rtol=1e-4, atol=1e-7)
    sol=solve_ivp(lambda x,y: odefuncexfull(x,y,r,DvYSZ, DvCGO, x0, KehYSZ, KehCGO,\
        bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO,epsilon,inindex,\
        idindex, F, i), [L, 0], y0, method=methodRK, t_eval=xspan,\
        dense_output=False, events=None, vectorized=False, rtol=comtol, atol=1e-6)
    

    #sol=odeint(lambda x,y: odefuncexfull(x,y,r,DvYSZ, DvCGO, x0, KehYSZ, KehCGO,\
     #   bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO,epsilon,inindex,\
     #   idindex, F, i),  y0, xspan)
        #sol=BDF(lambda x,y: odefuncexfull(x,y,r,DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO,epsilon,inindex, idindex, F, i), xspan[0], y0, xspan[-1])
    #[t,y,e]=rk45(lambda x,y: odefuncexfull(x,y,r,DvYSZ, DvCGO, x0, KehYSZ, \
    #        KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ,\
     #        DeCGO,conYSZ,conCGO,epsilon,inindex, idindex, F, i),\
    #             [xspan[0], xspan[-1]], np.array(y0), 80000)  
    Cexp=sol.y[0]
    #Cexp=sol[:,0]
    #Cexp=np.exp(lnCexp)
    xp=sol.t
    #xp=xspan
    #Cexp=y
    #xp=t
    #print('LengthCe=', len(Cexp))
    #print('Lengthxp=', len(xp))
    
    iv=i/(1+1/r)

 ###############solve Galvani potential \phi ###########
    ic = 0
#opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
#[xx,phi] = ode45(@(x,y) odephipf(r,DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO, conYSZ, conCGO,epsilon,xp,Cexp,x,y), xspan, ic);
#phi=phi-interp1(xx,phi,L/2);
    ####sol=solve_ivp(lambda x,y: odephipf(x,y,r,DvYSZ, DvCGO, x0, KehYSZ, KehCGO,\
    ####    bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO, conYSZ, conCGO, epsilon, xp, Cexp,\
    ####    idindex, inindex, R, T,  F, i, IndexInitial), [xspan[0], xspan[-1]], [0],\
    ####    method='RK45', t_eval=xspan, dense_output=False, events=None, \
    ####    vectorized=False, rtol=comtol, atol=1e-6)
    ####phi=sol.y[0]
    ####mt=sol.t
    ####print('lenghtxp=',len(xspan))
    ####print('Lengthfp=',len(mt))
    
    ####if IndexInitial==1:
    ####    reverse_xspan=mt[::-1]
    ####    reverse_phi=phi[::-1]
    ####    phimid=np.interp(L/2,reverse_xspan,reverse_phi)
    ####    phi=phi-phimid*np.ones(len(mt))
    ####elif IndexInitial==-1:
    ####    phimid=np.interp(L/2,mt,phi)
    ####    phi=phi-phimid*np.ones(len(mt))

    #pivolt=np.add(-(R*T/F)*np.log(Cexp),phi)
    #electrop=np.subtract(pivolt,phi)
    #muebar=np.subtract(R*T*np.log(Cexp),F*phi)
    cationYSZ=8/(AYSZ**3*10**(-24))/Avorgacon/conYSZ/(0.5*(bYSZ+1-KehYSZ/1))
    #print('cation',cationYSZ)
    cationGDC=8/(AGDC**3*10**(-24))/Avorgacon/conCGO/(0.5*(bCGO+1-KehCGO/1))
    #print('shapecation',type(cationYSZ))
    #print('shapexp',type(xp))
    #print('shapeCexp',type(Cexp))
    #print('sizeCexp',Cexp.ndim)
    #muhbar=np.zeros(len(xp))
    Chxp=np.zeros(len(xp))
    Cvxp=np.zeros(len(xp))
    Coxp=np.zeros(len(xp))
    #print('sizeCoxp',Coxp.ndim)
    for kk in range(0,len(xp),1):
     #   muhbar[kk]=R*T*np.log(IntepoAB(KehYSZ*conYSZ/conYSZh,KehCGO*conCGO/conCGOh,x0,epsilon,xp[kk])/Cexp[kk])+F*phi[kk]
        Chxp[kk]=IntepoAB(KehYSZ,KehCGO,x0,epsilon,xp[kk])/Cexp[kk]
        Cvxp[kk]=(IntepoAB(bYSZ,bCGO,x0,epsilon,xp[kk])+Cexp[kk]-\
                  IntepoAB(KehYSZ,KehCGO,x0,epsilon,xp[kk])/Cexp[kk])\
            /(IntepoAB(bYSZ,bCGO,x0,epsilon,xp[kk])+1\
              -IntepoAB(KehYSZ,KehCGO,x0,epsilon,xp[kk]))    
        Coxp[kk]=(IntepoAB(cationYSZ,cationGDC,x0,epsilon,xp[kk])-Cvxp[kk])\
            /(IntepoAB(cationYSZ,cationGDC,x0,epsilon,xp[kk])-1)

    if idindex==1 and inindex==1:
       if po2index==0:
           po2xp=np.exp(-(R*T/F)*np.log(Cexp)*4*F/R/T)
           #po2xp=Cexp**(-4)
       elif po2index==1:
           intarray=np.array(Cexp)
           po2xp=np.power((intarray.flatten()),(-4))*np.power((Cvxp),(-2))*np.power((Coxp),2)
           #po2xp=np.exp(-(R*T/F)*np.log(Cexp)*4*F/R/T)*np.power((np.power((Cvxp),(-1))*np.power((Coxp),1)),2)
    else:
      po2xp=np.exp(-(R*T/F)*np.log(Cexp)*4*F/R/T)

    xxspan=np.array(xp)/L
    err=0
    #plt.plot(xxspan,po2xp) 
    #plt.yscale("log")
    #plt.show
    return xxspan, po2xp, err



                          
        
