import math
def DiffIntepoAB(A,B,x0,epsilon,x):

#y=B*0.5*(1+tanh((x-x0)/(sqrt(2)*epsilon)))+A*0.5*(1+tanh((x0-x)/(sqrt(2)*epsilon)));

#y=B*0.5*(1+tanh((x-x0)/(sqrt(2)*epsilon)))+A;
    y=-(math.sqrt(2)*(B - A)*(math.tanh((math.sqrt(2)*(x - x0))/(2*epsilon))**2 - 1))/(4*epsilon);

    return y
