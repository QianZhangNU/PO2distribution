import numpy as np

def OxygenParT(T,PH2,PH2O):

    y=1.295-0.243*((T+273.15)-300)/700
    PO2=0.2/np.exp(y*4*96400/8.314/(T+273.15))
    PO2=PO2*(PH2O/PH2)**2
    
    return PO2
