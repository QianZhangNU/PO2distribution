# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 15:27:01 2020

@author: qianj
"""
import numpy as np

class LeftLayer:
    def __init__(self, T_temperature):
        self.T_temperature=T_temperature
        T=T_temperature+273.15
        Avorgacon=6.02*1e+23
        eVoverkT=1.0/8.618*1e+5/T
        R=8.314  #J/(mol*K)
        F=96485
        sigmaion=163*np.exp(-0.79*eVoverkT)
        #self.sigmaion=163*np.exp(-0.79*eVoverkT)
        self.De = 230.0*np.exp(-2*eVoverkT)
        self.Dh = 0.23*np.exp(-1.15*eVoverkT)
        Keh=(1.04*10**23*np.exp(-1.99*eVoverkT))/Avorgacon*1.72*10**21*np.exp(-0.62*eVoverkT)/Avorgacon
        con=(1.04*10**23*np.exp(-1.99*eVoverkT))/Avorgacon #mol/cm^3
        conh=1.72*10**21*np.exp(-0.62*eVoverkT)/Avorgacon
        self.Keh=Keh/con/con
        self.A=5.1385
        AYSZ=5.1385
        CY=2*8/(AYSZ**3*10**(-24))/Avorgacon*0.08/2.16
    #CY=0.0078288; %mol cm-3
        self.Dv=sigmaion*R*T/4/F/F/(CY/2)
    #pause
        cGd=CY #mol/cm3;
        cGd=cGd/con
        self.b=cGd
        self.con=con


class RightLayer:
    def __init__(self,T_temperature):
        T=T_temperature+273.15
        Avorgacon=6.02*1e+23
        eVoverkT=1.0/8.618*1e+5/T
        R=8.314  #J/(mol*K)
        F=96485
        muh0=2.4*1e+6*np.exp(-1.45*1.0/8.618*1e+5/1073.15)/1073.15/(1.602*1e-19)/(9.26e+18)/np.exp(-1.08*1.0/8.618*1e+5/1073.15)*1073.15
        muhCGO=muh0/T*np.exp(-1.08*eVoverkT)
        conCGO=4*1e+9*np.exp(-2.51*eVoverkT)/(1.602*1e-19)/(10**(-4301/T+4.1947))/Avorgacon #mol/cm^3
        conCGOh=2.4*1e+6*np.exp(-1.45*eVoverkT)/T/muhCGO/(1.602*1e-19)/Avorgacon
        conCGOh=2.4*1e+6*np.exp(-1.45*eVoverkT)/(1392*np.exp(-1.08*eVoverkT))/(1.602*1e-19)/Avorgacon
        DvCGO=10**(-3405/T-1.8995)
        DeCGO=10**(-4301/T+0.1301) #cm^2/s
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
        self.De=DeCGO
        self.Dh=DhCGO
        self.Dv=DvCGO
        self.Keh=KehCGO
        self.A=AGDC
        self.b=cGdCGO
        self.con=conCGO

#YSZ=LeftLayer(800)
#print('DeYSZ=',YSZ.De)
    