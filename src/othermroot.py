# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 14:57:58 2020

@author: qianj
"""


import numpy as np
from fzeroforrf import fzeroforrf

def othermroot(DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, rinitial, valrtol, valatol):
    #xx=np.arange(xleft,xright+(xright-xleft)/N,(xright-xleft)/N) 
    #fxx=[fzeroforrf(x, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i) for x in xx]
    #fxx=np.array(fxx)
    #cs=np.multiply(fxx,np.roll(fxx,-1))
    #xc=[]
     #   for kk in range(0,len(cs)):
     #       if cs[kk]<=0:
     #           xc.append(xx[kk])
                
     #   xc = np.array(xc)
     xleft=rinitial
     xright=rinitial
     for n in range(0,800,1):
         interval=np.abs(2.0*rinitial)/1600
         xleft=xleft-interval
         xright=xright+interval
         fleft=fzeroforrf(xleft, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol)
         fright=fzeroforrf(xright, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol)
         if (fleft*fright)<=0:
             a=xleft
             b=xright
         else:
             continue
         
     return a, b
         