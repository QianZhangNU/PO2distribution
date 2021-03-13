 import ParaCond as PC

class InitialGuess:
    def __init__(self,c0,cL,i_current,L,x_portion,T_temperature,RpH,RpO,PO2EH,PO2EO):
        self.c0=c0
        self.cL=cL
        self.i_current=i_current
        self.L=L
        self.x_portion=x_portion
        self.T_temperature=T_temperature
        self.RpH=RpH

    def rinitwguess(self,formerindex):      
        Estfactor=factorfront(i_current,L,x_portion,T_temperature)
    po2oeled=np.exp(np.log(PO2EO)-2*R*T/(2*F)*np.arcsinh(i/(2*R*T/2/F/RpO))*4*F/R/T)
    po2guess=po2oeled*np.exp(Estfactor)
    cLguess=np.power(po2guess,(-1/4))
    #print('testinf=',np.sinh(i*RpO*F/R/T))
    print('Est=',Estfactor)
    print('po2eled=',po2oeled)
    print('po2guess=',po2guess)
    print('CL=',cL)
    print('clguess=',cLguess)
    ##### Solve r and Ce ################
    print("test IndexInitial = " + str(IndexInitial))
    if formerindex==1:
       #rinitial=potenfindr(i,RpH,RpO,PO2EH,PO2EO,L*(1-x_portion),T-273.15,1,c0,cLguess)
       rinitial=potenfindr(i,RpH,RpO,PO2EH,PO2EO,L,T-273.15,1,c0,cL)
        #print('rinitial1=', rinitial)
       #rinitial=otherzero(DvYSZ, DvCGO, x0, KehYSZ, KehCGO, bYSZ, bCGO, DhYSZ, DhCGO, DeYSZ, DeCGO,conYSZ,conCGO, epsilon, IndexInitial, c0, cL, L, inindex, idindex, F, i, rinitial, 2e+4)
       #print(rinitial)
    elif formerindex==-1:
       rinitial=-0.972
    else:
       rinitial=rformer





