import numpy as np
def PO2WCV(x,KehYSZ,bYSZ,conYSZ,A,PO2):
    Avorgacon=6.02*1e+23
    Cv=(bYSZ+x-KehYSZ/x)/(bYSZ+1-KehYSZ/1)
    Canoin=8/(A**3*10**(-24))/Avorgacon/conYSZ/(0.5*(bYSZ+1-KehYSZ/1))
    Co=(Canoin-Cv)/(Canoin-1)
    Ccation=4/(A**3*10**(-24))/Avorgacon/conYSZ
  #  Czr=(Ccation-x-bYSZ)/(Ccation-1-bYSZ);

    y=(-2*np.log((x)**2)-np.log((Cv)**2)+np.log((Co)**2))-np.log(PO2)
    return y
