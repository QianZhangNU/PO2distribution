# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 09:22:56 2020

@author: qianj
"""
import math
from fzeroforrf import fzeroforrf

def zerointerval(DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol, xinitial,methodRK):
     xleft=xinitial
     xright=xinitial
      
     if xinitial != 0:
         dx=0.0283*xinitial
     else:
         dx=0.0283
         
     xleft=xleft-dx
     xright=xright+dx
     fleft=fzeroforrf(xleft, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol,methodRK)
     fright=fzeroforrf(xright, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol,methodRK)
     if fleft*fright<=0:
         a=xleft
         b=xright
     
     
     dx=(math.sqrt(2)-1)*dx
    
     xleft=xleft-dx
     xright=xright+dx
     fleft=fzeroforrf(xleft, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol,methodRK)
     fright=fzeroforrf(xright, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol,methodRK)
     if fleft*fright<=0:
         a=xleft
         b=xright
     
     mm=1;
     while (fleft*fright>0) and (xleft>0) and (xright<2*xinitial) and (mm<50):
         dx=math.sqrt(2)*dx
         mm=mm+1
         xleft=xleft-dx
         xright=xright+dx
         fleft=fzeroforrf(xleft, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol,methodRK)
         fright=fzeroforrf(xright, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, valrtol, valatol,methodRK)
         if fleft*fright<=0:
             break
         elif abs(fleft+1)<=1e-16: 
             if abs(fright+1)>1e-16:
                 break
             else: 
                 continue
         elif abs(fright+1)<=1e-16:
             if abs(fleft+1)>1e-16:
                 break
             else:
                 continue
         #print('fleft=',fleft)
         #print('fright=',fright)
         
     a=xleft
     b=xright
     if abs(fleft)>=abs(fright):
         c=xright
     else:
         c=xleft
         
     return a, b, c, mm
            
    