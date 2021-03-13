import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import BDF
from odefuncexfull import odefuncexfull
from RK import r8_rkf45
from rk45 import rk45
#from assimulo.solvers import CVode
#from assimulo.solvers import LSODAR
#from assimulo.problem import Explicit_Problem
def fzeroforrf(tt,DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon,IndexInitial, c0, cL, L,inindex, idindex, F, i, valrtol, valatol,methodRK):
   # global c0, cL, L

    if IndexInitial==1:
       y0=cL
       #xspan=np.arange(L,-epsilon/4,-epsilon/4)
       
       xspan=L-np.linspace(0,L,num=8e+4+1)
       #print('xspan[0]',xspan[0])
       #print('xspan[len]',xspan[len(xspan)-1])
       sol=solve_ivp(lambda x,y: odefuncexfull(x,y,tt,DvYSZ, DvCGO, x0, KehYSZ, \
            KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ,\
             DeCGO,conYSZ,conCGO,epsilon,inindex, idindex, F, i),\
                  [L, 0], y0, method=methodRK, t_eval=xspan,\
                      dense_output=False, rtol=valrtol, atol=valatol)
       #####model = Explicit_Problem(lambda x,y: odefuncexfull(x,y,tt,DvYSZ, DvCGO, x0, KehYSZ, KehCGO,\
       #####bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO,epsilon,inindex,\
       #####idindex, F, i), y0, L)
       #####sim = CVode(model)
       #####sim.rtol=1.e-3
       #####sim.atol=1.e-6
       #####tfinal = 0.0        #Specify the final time
       #####sim.backward = True
       #####xp, yp = sim.simulate(tfinal)
       # y, yp, t, flag  = r8_rkf45(f, 1, y, yp, t, tout, relerr, abserr, flag)
       #sol=solve_ivp(lambda x,y: odefuncexfull(x,y,tt,DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO,epsilon,inindex, idindex, F, i), [xspan[0], xspan[-1]], y0, method='RK45', t_eval=xspan, dense_output=False, events=None, vectorized=False)
       #sol=BDF(lambda x,y: odefuncexfull(x,y,tt,DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO,epsilon,inindex, idindex, F, i), L, y0, 0)
       y=sol.y[0]
       #[t,y,e]=rk45(lambda x,y: odefuncexfull(x,y,tt,DvYSZ, DvCGO, x0, KehYSZ, \
       #     KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ,\
        #     DeCGO,conYSZ,conCGO,epsilon,inindex, idindex, F, i),\
         #        [xspan[0], xspan[-1]], np.array(y0), 80000)   
                    
 #  [x,y]=ode45(@(x,y) odefuncexfull(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon,x,y), xspan, y0);
       #xp, yp = sim.simulate(tfinal)
       yy=(y[len(y)-1]-(c0))/c0
       #####yy=(y[len(y)-1]-(c0))/(c0)
    elif IndexInitial==-1:
         y0=c0
         xspan=np.arange(0, L+epsilon/4, epsilon/4)
#   [x,y]=ode45(@(x,y) odefuncexfull(tt, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon,x,y), xspan, y0);
         sol=solve_ivp(lambda x,y: odefuncexfull(x,y,tt,DvYSZ, DvCGO, x0, KehYSZ, \
              KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,\
               conYSZ,conCGO,epsilon,inindex, idindex, F, i),\
                [xspan[0], xspan[-1]], y0, method='RK23', \
                  dense_output=False,  rtol=valrtol, atol=volatol)
         #sol=solve_ivp(lambda x,y: odefuncexfull(x,y,tt,DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO,epsilon,inindex, idindex, F, i), [xspan[0], xspan[-1]], y0, method='RK45', t_eval=xspan, dense_output=False, events=None, vectorized=False)
         #sol=BDF(lambda x,y: odefuncexfull(x,y,tt,DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO,epsilon,inindex, idindex, F, i), 0, y0, L)       
         y=sol.y[0]
         yy=(y[len(y)-1]-(cL))/(cL)


    return yy
