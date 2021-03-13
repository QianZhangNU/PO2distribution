import numpy as np
from IntepoAB import IntepoAB
from DiffIntepoAB import DiffIntepoAB

def odephipf(x,y,r,DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO, conYSZ, conCGO, epsilon,xx,cexx,idindex, inindex, R, T,  F, i, IndexInitial):

  #  global  idindex, inindex, R, T,  F, i 

    if IndexInitial==1:
        reverse_xx=xx[::-1]
        reverse_cexx=cexx[::-1]
        Ce=np.interp(x,reverse_xx,reverse_cexx)
  #      phi=phi-phimid*np.ones(len(xspan))
    elif IndexInitial==-1:
        Ce=np.interp(x,xx,cexx)
        

#Ce=interp1(xx,cexx,x)
#dydx=-iv*(Ce.*R*T*(Dv*Keh*id + Dh*Keh*r + Ce.^2*Dv*id + Ce.^2*De*r))/(Dv*F^2*r*(De*Ce.^2 + Dh*Keh)*(Keh*id - 2*Keh*in + Ce.^2*id + 2*Ce.^2*in + 2*CGd.*Ce))

#dydx=-iv*(Ce.*R*T*(IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*id + IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*r + Ce.^2*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*id + Ce.^2*IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*r))/(IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F^2*r*(IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*Ce.^2 + IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x))*(IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*id - 2*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*in + Ce.^2*id + 2*Ce.^2*in + 2*IntepoAB(bYSZ,bCGO,x0,epsilon,x).*Ce));

    dydx=i*(IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*idindex +IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*r + Ce**2*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*idindex + Ce**2*IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*r)

    dydx=dydx-Ce*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*DiffIntepoAB(KehYSZ*conYSZ,KehCGO*conCGO,x0,epsilon,x)*idindex*(IntepoAB(DeYSZ,DeCGO,x0,epsilon,x) - IntepoAB(DhYSZ,DhCGO,x0,epsilon,x))*(r + 1)

    dydx=dydx+IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*DiffIntepoAB(bYSZ*conYSZ,bCGO*conCGO,x0,epsilon,x)*idindex*(r+1)*(IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*Ce**2+IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x))

    dydx=-Ce*R*T*dydx

    dydx=dydx-IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*(r+1)*Ce**2*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*R*T*DiffIntepoAB(conYSZ,conCGO,x0,epsilon,x)*idindex*(IntepoAB(DeYSZ,DeCGO,x0,epsilon,x) -IntepoAB(DhYSZ,DhCGO,x0,epsilon,x))

    dydx=dydx/(IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F**2*IntepoAB(conYSZ,conCGO,x0,epsilon,x)*(r+1)*(IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*Ce**2 + IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x))*(IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*idindex- 2*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*inindex + Ce**2*idindex + 2*Ce**2*inindex + 2*IntepoAB(bYSZ,bCGO,x0,epsilon,x)*Ce))

    return dydx
