import math

def IntepoAB(A,B,x0,epsilon,x):

    y=B*0.5*(1+math.tanh((x-x0)/(math.sqrt(2)*epsilon)))+A*0.5*(1+math.tanh((x0-x)/(math.sqrt(2)*epsilon)))

#y=B*0.5*(1+tanh((x-x0)/(sqrt(2)*epsilon)))+A;


    return y
