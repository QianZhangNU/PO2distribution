import numpy as np
import math

def bdcsym(RpH,RpO,PO2EH,PO2EO,i,R,T,F):
    PO2H=np.exp(np.log(PO2EH)+2*R*T/(2*F)*np.arcsinh(i/(2*R*T/2/F/RpH))*4*F/R/T);
    PO2O=np.exp(np.log(PO2EO)-2*R*T/(2*F)*np.arcsinh(i/(2*R*T/2/F/RpO))*4*F/R/T);


    c0=np.power(PO2H,(-1/4))
    cL=np.power(PO2O,(-1/4))
     
    return c0, cL
