import numpy as np

def bdcWCV(RpH,RpO,EOCVH,EOCVO,i,R,T,F):
#diffvH=diffveqH+i*RctH;
#diffvO=diffveqO-i*RctO;
#global R T F
    PO2H=np.exp(EOCVH*F/R/T+2*R*T/(2*F)*np.arcsinh(i*RpH/(2*R*T/2/F))*F/R/T)
    PO2O=np.exp(EOCVO*F/R/T-2*R*T/(2*F)*np.arcsinh(i*RpO/(2*R*T/2/F))*F/R/T)


    c0=PO2H**-1
    cL=PO2O**-1

    return c0, cL
