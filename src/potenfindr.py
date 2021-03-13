import numpy as np
import math
from scipy.optimize import root
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
from initialguessboverc3 import intialguessboverc3
from findBoverC import findBoverC

def  potenfindr(i_current,RpH,RpO,PO2EH,PO2EO,L_length,T_temperature,Indicator,c0,cL):
  #   global R, i, F, T,  sigmaion, De, A, cL, c0, L, eVoverkT, con,  B, BoverC 

#parameters
     R=8.314 #J/(mol*k)
     T=T_temperature+273.15
     i=i_current
     F=96485.0 #C/mol
     Avorgacon=6.02*1e+23 #mol^-1
#RpH=Rp #OHM*cm^2
#RpO=Rp #OHM*cm^2
     L=L_length
     eVoverkT=1.0/8.618*1e+5/T
     sigmaion=163*np.exp(-0.79*eVoverkT) #OHM^-1cm^-1
     De=230.0*np.exp(-2*eVoverkT) #cm^2/s
     Dh=0.23*np.exp(-1.15*eVoverkT)
     Keh=(1.04*1e+23*np.exp(-1.99*eVoverkT))/Avorgacon*1.72*1e+21*np.exp(-0.62*eVoverkT)/Avorgacon
     con=(1.04*1e+23*np.exp(-1.99*eVoverkT))/Avorgacon #mol/cm^3
     A=Dh*Keh/De/con/con


#boundary conditions of concentration
# if (Indicator==1)
#    [c0, cL]=bdcsym(RpH,RpO,PO2EH,PO2EO,i,R,T,F)
# elseif (Indicator==-1)
#    [c0, cL]=bdcsymH2(RpH,RpO,0.97,0.97,0.03,0.03,i,R,T,F)
# else
#    disp('Indicator has to be either 1 or -1.');
#    pause
# end
# cL=cL*conCGO/conYSZ;
#pause
     rangeofboverc=1/math.sqrt(4*A)
      #fun=@findBoverC;

#set initial guess
     y0=intialguessboverc3(rangeofboverc,2e+5, R, T, F, i, sigmaion, De, A, cL, c0, L, con)
    # print(len(y0))
#pause
     for ii in range(len(y0)-1,0,-1):
#for ii=1:1:length(y0)
         x0=y0[ii]
         [xx, infodict, exitflag, mesg] = fsolve(lambda x: findBoverC(x, R, T, F, i, sigmaion, De, A, cL, c0, L, con), x0, fprime=None, full_output=1)
         if exitflag==1:
             break
         else:
             continue

     BoverC=xx
  #   print(xx)
# B=1/(-sigmaion*R*T/F/i*(1-F*F*De*con/sigmaion/R/T/BoverC));
# iion=-sigmaion*R*T*B/F;
# iionoieh=iion/(i-iion);
     iionoieh=BoverC*(-R*T*sigmaion/F/F/De/con)
   #  print(iionoieh)
     rinitial=np.asscalar(iionoieh)
   #  print(type(rinitial))
     
     return rinitial


        






