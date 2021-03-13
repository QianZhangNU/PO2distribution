# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 18:51:18 2020

@author: qianj
"""
import numpy as np

def factorfront(i_current,L,x_portion,T_temperature):

   R=8.314  #J/(mol*k)
   F=96485.0 #C/mol
   Avorgacon=6.02*1e+23 #mol^-1


#%%%%%%%Operating Conditions of the Cell%%%%%%%%%%%%%

   T=T_temperature+273 #K
   i=i_current #Acm^-2
   eVoverkT=1.0/8.618*1e+5/T

#%%%%% PARAMETERS IN YSZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%% PARAMETERS IN YSZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%

   sigmaionYSZ=163*np.exp(-0.79*eVoverkT)
   DeYSZ=230.0*np.exp(-2*eVoverkT)
   DhYSZ=0.23*np.exp(-1.15*eVoverkT)
   KehYSZ=(1.04*10**23*np.exp(-1.99*eVoverkT))/Avorgacon*1.72*10**21*np.exp(-0.62*eVoverkT)/Avorgacon
   conYSZ=(1.04*10**23*np.exp(-1.99*eVoverkT))/Avorgacon
   conYSZh=1.72*10**21*np.exp(-0.62*eVoverkT)/Avorgacon
   KehYSZ=KehYSZ/conYSZ/conYSZ
   AYSZ=5.1385
   CY=2*8/(AYSZ**3*10**(-24))/Avorgacon*0.08/2.16
#CY=0.0078288; %mol cm-3
   DvYSZ=sigmaionYSZ*R*T/4/F/F/(CY/2)*1
#pause
   cGdYSZ=CY
   cGdYSZ=cGdYSZ/conYSZ
   bYSZ=cGdYSZ
   sigmaChYSZ=DhYSZ*F**2*conYSZh/R/T


#%%%%%%%PARAMETERS IN CGO%%%%%%%%%%%%%%%%%%%%%%%%%%%
   muh0=2.4*1e+6*np.exp(-1.45*1.0/8.618*10**5/1073.15)/1073.15/(1.602*1e-19)/(9.26e+18)/np.exp(-1.08*1.0/8.618*10**5/1073.15)*1073.15
   muhCGO=muh0/T*np.exp(-1.08*eVoverkT)
   conCGO=4*1e+9*np.exp(-2.51*eVoverkT)/(1.602*1e-19)/(10**(-4301/T+4.1947))/Avorgacon
   conCGOh=2.4*1e+6*np.exp(-1.45*eVoverkT)/T/muhCGO/(1.602*1e-19)/Avorgacon
#%conCGOh=2.4*1e+6*exp(-1.45*eVoverkT)/(1392*exp(-1.08*eVoverkT))/(1.602*1e-19)/Avorgacon*factortorR;
   DvCGO=10**(-3405/T-1.8995)*1
   DeCGO=10**(-4301/T+0.1301)
   DhCGO=muh0*R/F*np.exp(-1.08*eVoverkT)
#%DhCGO=0.12*exp(-1.08*eVoverkT)*factortorR;
   KehCGO=conCGO*conCGOh
   KehCGO=KehCGO/conCGO/conCGO
#%pause
   AGDC=5.418
   cGdCGO=2*8/(AGDC**3*10**(-24))/Avorgacon*0.05/2.0
   sigmaionCGO=DvCGO*4*F*F*(cGdCGO/2)/R/T
#%cGdCGO=0.0041789*factortorR; %mol/cm3;
#%conCGO-conCGOh
   cGdCGO=cGdCGO/conCGO
   bCGO=cGdCGO
   sigmaChCGO=DhCGO*F**2*conCGOh/R/T

    #leftOPO2=log10(po2);
    #rightOPO2=leftOPO2;
    #[PHOCV,c0I,cce0,CvGDC,EA]=PO2VSCEGDC(T_temperature,leftOPO2,rightOPO2);
    #[PHOCV,c0I,cce0,CvYSZ,EA]=PO2VSCE(T_temperature,leftOPO2,rightOPO2);

   factorlnpoi=i*L*x_portion*4*F/R/T*(sigmaChYSZ/(sigmaionYSZ)-sigmaChCGO/(sigmaionCGO))*(1/sigmaChCGO)
    
   return factorlnpoi