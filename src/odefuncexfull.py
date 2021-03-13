import numpy as np
import math
from IntepoAB import IntepoAB
from DiffIntepoAB import DiffIntepoAB

def odefuncexfull(x,y,r,DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO,epsilon,inindex, idindex, F, i):

  #  global inindex, idindex, F, i 
#DvYSZ DvCGO x0 KehYSZ KehCGO bYSZ bCGO DhYSZ DhCGO DeYSZ DeCGO 
# dydx=(2*y^3*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*i*in - y^3*IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*i*r + 2*IntepoAB(bYSZ,bCGO,x0,epsilon,x)*y^2*Dv*ii - 2*y*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*i*in - y*IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*i*r);
# dydx=y*dydx+y*con*( 2*IntepoAB(bYSZ,bCGO,x0,epsilon,x)*y*IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*DiffIntepoAB(KehYSZ,KehCGO,x0,epsilon,x) + IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*DiffIntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*id - 2*IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*DiffIntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*in + y^2*IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*DiffIntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*id + 2*y^2*IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*DiffIntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*in + y^2*IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*DiffIntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*id*r + 2*y^2*IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*DiffIntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*in*r + 2*IntepoAB(bYSZ,bCGO,x0,epsilon,x)*y*IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*DiffIntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*r + IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*DiffIntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*id*r - 2*IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*DiffIntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*in*r);
# dydx=dydx/(IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*(IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*y^2 + IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x))*(r + 1)*(IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*id - 2*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*in + y^2*id + 2*y^2*in + 2*IntepoAB(bYSZ,bCGO,x0,epsilon,x)*y));
# dydx=dydx/con;
    #print('yCe=', y)
    #y=np.exp(y)

    dydx=-y*i*(2*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*inindex - 2*IntepoAB(bYSZ,bCGO,x0,epsilon,x)*y*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x) + IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*r - 2*y**2*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*inindex + y**2*IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*r)

    dydx=dydx+IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*DiffIntepoAB(KehYSZ*conYSZ,KehCGO*conCGO,x0,epsilon,x)*(r + 1)*(2*IntepoAB(bYSZ,bCGO,x0,epsilon,x)*y*IntepoAB(DhYSZ,DhCGO,x0,epsilon,x) +IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*idindex - 2*IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*inindex + y**2*IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*idindex + 2*y**2*IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*inindex)

    dydx=dydx-y**2*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*DiffIntepoAB(conYSZ,conCGO,x0,epsilon,x)*(r + 1)*(2*IntepoAB(bYSZ,bCGO,x0,epsilon,x)*y*IntepoAB(DeYSZ,DeCGO,x0,epsilon,x) + IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*idindex - 2*IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*inindex + y**2*IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*idindex + 2*y**2*IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*inindex)

    dydx=dydx-y*IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*DiffIntepoAB(bYSZ*conYSZ,bCGO*conCGO,x0,epsilon,x)*(r + 1)*idindex*(IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*y**2+IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x))

    dydx=y*dydx

    dydx=dydx/(IntepoAB(DvYSZ,DvCGO,x0,epsilon,x)*F*(IntepoAB(DeYSZ,DeCGO,x0,epsilon,x)*y**2 + IntepoAB(DhYSZ,DhCGO,x0,epsilon,x)*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x))*(r + 1)*(IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*idindex - 2*IntepoAB(KehYSZ,KehCGO,x0,epsilon,x)*inindex + y**2*idindex + 2*y**2*inindex + 2*IntepoAB(bYSZ,bCGO,x0,epsilon,x)*y))

    dydx=dydx/IntepoAB(conYSZ,conCGO,x0,epsilon,x)
    
    #dydx=y*dydx

    return dydx
