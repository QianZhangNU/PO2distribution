import cmath
import math
import numpy as np


def findBoverC(x,R, T, F, i, sigmaion, De, A, cL, c0, L, con):
 #   global R, T, F, i, sigmaion, De, A, cL, c0, L, con
#x=x/1e+16;
#x=x/1e+8;
#x=x/1e+14;
#x=x*-1.062151850602328e-13;
#ff=log((2*cL*x+1+sqrt(1-4*A*x*x))/(2*c0*x+1+sqrt(1-4*A*x*x)));
#ff=ff-log((2*cL*x+1-sqrt(1-4*A*x*x))/(2*c0*x+1-sqrt(1-4*A*x*x)));
    if (1-4*A*x*x)>0:
        if x==0:
           ff=1
        else:
           fz=4*(c0*cL+A)*x*(x)+2*x*(c0+cL)+2*x*(c0-cL)*math.sqrt(1-4*A*(x)*(x))
           fx=4*(c0*cL+A)*x*(x)+2*x*(c0+cL)+2*x*(cL-c0)*math.sqrt(1-4*A*(x)*(x))
#fx=double(fx);
           ff=fz/(fx+1e-60)
#  (2*cL*x+1+sqrt(1-4*A*x*x))
#  (2*c0*x+1-sqrt(1-4*A*x*x))
#  (2*cL*x+1-sqrt(1-4*A*x*x))
#  (2*c0*x+1+sqrt(1-4*A*x*x))
#end

        ff=np.log(abs(ff))
    
#y=-L-sigmaion*R*T/F/i*(1-F*F*De*con/sigmaion/R/T/x)*(log(cL/c0)+1/sqrt(1-4*A*x*x)*ff);
#y
    if (1-4*A*x*x)>0:
#y=-L-sigmaion*R*T/F/i*(1-F*F*De*con/sigmaion/R/T/x)*(log(cL/c0)+1/sqrt(1-4*A*x*x)*ff);
        y=-L-sigmaion*R*T/F/i*(1-F*F*De*con/sigmaion/R/T/x)*(np.log(cL/c0)+1/math.sqrt(1-4*A*x*x)*ff)
    elif (1-4*A*x*x)<0:
#y=-L-sigmaion*R*T/F/i*(1-F*F*De*con/sigmaion/R/T/x)*(log(cL/c0)-2/sqrt(4*A*x*x-1)*(atan((2*x*cL+1)/(sqrt(4*A*x*x-1)))-atan((2*x*c0+1)/(sqrt(4*A*x*x-1)))));
        y=-L-sigmaion*R*T/F/i*(1-F*F*De*con/sigmaion/R/T/x)*(np.log(cL/c0)-2/math.sqrt(4*A*x*x-1)*(math.atan((2*x*cL+1)/(math.sqrt(4*A*x*x-1)))-math.atan((2*x*c0+1)/(math.sqrt(4*A*x*x-1)))))
    else:
#y=-L-sigmaion*R*T/F/i*(1-F*F*De*con/sigmaion/R/T/x)*(log(cL/c0)+1/x*(1/(cL+0.5/x)-1/(c0+0.5/x)));
        y=-L-sigmaion*R*T/F/i*(1-F*F*De*con/sigmaion/R/T/x)*(np.log(cL/c0)+(1/(x*cL+0.5)-1/(x*c0+0.5)))
#end

#y=y^(5);
    return y


