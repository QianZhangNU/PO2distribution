# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:59:59 2020

@author: qianj
"""


import numpy as np
import math
from fzeroforrf import fzeroforrf

def otherzero(DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, rinitial, N1):
    N1=N1+1
    #A=DhYSZ/DeYSZ*KehYSZ
    #rangeofboverc=1/math.sqrt(4*A)
    #xleft=-rangeofboverc
    #xright=rangeofboverc
    xleft=0
    xright=2*rinitial
    for n in range(0,10,1):
        xx=np.arange(xleft,xright+(xright-xleft)/N1,(xright-xleft)/N1)
        fxx=[]
        fxx=[fzeroforrf(x, DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i) for x in xx]
    #fxx=np.array(fxx)
        cs=np.multiply(fxx,np.roll(fxx,-1))       
        xc=[]
        for kk in range(0,len(cs)):
            if cs[kk]<=0:
                xc.append(xx[kk])
                
        xc = np.array(xc)
#xc = xxx(cs <= 0);   
        if np.all(xc==0):
      #     k=0
           xxleft=[]
           xxright=[]
           m=[]
           for ii in range(1,len(xx)-1,1):
        #    secondd(ii)=(BB(ii)-BB(ii-1));
               if (fxx[ii]-fxx[ii-1])*(fxx[ii+1]-fxx[ii])<0:
                   xxleft.append(xx[ii-1])
                   xxright.append(xx[ii+1])
                   m.append(fxx[ii])
     #    break;
         #          k=k+1;
     #   end
   # end
           
           xxleft=np.array(xxleft)
           xxright=np.array(xxright)
           m=np.array(m)
   #    [~,I]=min(abs(m));
           I=np.argmin(abs(m)) 
           xleft=xxleft[I]
           print(xleft)
           xright=xxright[I]
           print(xright)
        else:
            y=xc
            print('lengthY=',len(y))
            break


    return y[1]