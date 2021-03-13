import numpy as np
from scipy.optimize import fsolve
from PO2WCV import PO2WCV
def EOCVWCV(KehYSZ,bYSZ,conYSZ,A,PO2,cI,R,T,F):
#fun=@(x) PO2WCV(x,KehYSZ,bYSZ,conYSZ,A,PO2);
#options = optimset('Display','iter'); 
#[xx fval exitflag output] = fzero(fun,cI,options);
  #  xx, infodict, ier, mesg=fsolve(PO2WCV,cI,(KehYSZ,bYSZ,conYSZ,A,PO2))
    xx=fsolve(PO2WCV,cI,(KehYSZ,bYSZ,conYSZ,A,PO2))
    EOCV=-R*T/F*np.log(xx)

    return EOCV
