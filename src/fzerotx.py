# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:18:13 2020

@author: Qian Zhang

"""
import numpy as np

def fzerotx(f,a,b,eps,*arg):
    #a=ab[1]
    #b=ab[2]
    fa=f(a,*arg)
    fb=f(b,*arg)
    c=a
    fc=fa
    d=b-c
    e=d
    while fb != 0:
        if fa*fb>0:
            a=c
            fa=fc
            d=b-c
            e=d
            
        if abs(fa)<abs(fb):
            c=b
            b=a
            a=c
            fc=fb
            fb=fa
            fa=fc
            
        m=0.5*(a-b)
        tol= 2.0*eps*max(abs(b),1.0)
        if (abs(m) <= tol) or (fb == 0.0):
           break
       
        if (abs(e) < tol) or (abs(fc) <= abs(fb)):
           d = m
           e = m
        else:
      # Interpolation
           s = fb/fc         
           if a == c:
         # Linear interpolation (secant)
              p = 2.0*m*s
              q = 1.0 - s
           else:
         # Inverse quadratic interpolation
              q = fc/fa
              r = fb/fa
              p = s*(2.0*m*q*(q - r) - (b - c)*(r - 1.0))
              q = (q - 1.0)*(r - 1.0)*(s - 1.0)
      
           if p > 0:
               q = -q
           else:
               p = -p
           
      # Is interpolated point acceptable
           if (2.0*p < 3.0*m*q - abs(tol*q)) and (p < abs(0.5*e*q)):
              e = d
              d = p/q
           else:
              d = m
              e = m
        
        # Next point
        c = b
        fc = fb
        if abs(d) > tol:
          b = b + d
        else:
          b = b - np.sign(b-a)*tol
   
        fb = f(b,*arg)
        
    return b
      
      
        
