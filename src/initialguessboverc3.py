import numpy as np
from findBoverC import findBoverC

def intialguessboverc3(rangeofboverc,N1,R, T, F, i, sigmaion, De, A, cL, c0, L, con):
#rangeofboverc
    N1=N1+1
    xleft=-rangeofboverc
    xright=rangeofboverc
    for n in range(0,10,1):
  #  n
#xx=-(rangeofboverc)/(10^(n-1)):2*rangeofboverc/(10^(n-1))/N1:(rangeofboverc)/(10^(n-1));]
      #  print('xleft = ' + str(xleft) + " , x_right = " + str(xright) + " , N1 = " + str(N1))  
        xx=np.arange(xleft,xright+(xright-xleft)/N1,(xright-xleft)/N1)
      #  print(len(xx))
        np.reshape(xx, (1,len(xx)))
        fxx=[findBoverC(x, R, T, F, i, sigmaion, De, A, cL, c0, L, con) for x in xx]
#fxx=arrayfun(@findBoverC,xx);
  #      ii=0
        fxx=np.array(fxx)
        xxx=[]
        fxxx=[]
  #      m=[]
     #   xxx = np.array([0 for _ in range(len(xx))])
     #   fxxx = np.array([0 for _ in range(len(xx))])
     #   m= np.array([0 for _ in range(len(xx))])
        xxleft = np.array([0 for _ in range(len(xx))])
        xxright = np.array([0 for _ in range(len(xx))])
        for iii in range(0,len(xx),1):
            if fxx[iii].imag==0:
                xxx.append(xx[iii])
                fxxx.append(fxx[iii])
      #          ii=ii+1 
# plot(xxx,fxxx,'*');
#pause
#hold on
#length(fxx)
#length(circshift(fxx,-1,1))
#cs = fxxx.*circshift(fxxx,-1,2);        % Product Negative At Zero-Crossings
        xxx=np.array(xxx)
        fxxx=np.array(fxxx)
        cs=np.multiply(fxxx,np.roll(fxxx,-1))
#plot(xxx,cs,'*');
#pause 
    #    cs.flatten()
    #    print(xxx.shape)
    #    print(xx.shape)
    #    print(cs.shape)
    #    print(fxxx.shape)
    #    xxx.flatten()
        xc=[]
        for kk in range(0,len(cs)):
            if cs[kk]<=0:
                xc.append(xxx[kk])
                
        xc = np.array(xc)
#xc = xxx(cs <= 0);   
        if np.all(xc==0):
      #     k=0
           xxleft=[]
           xxright=[]
           m=[]
           for ii in range(1,len(xxx)-1,1):
        #    secondd(ii)=(BB(ii)-BB(ii-1));
               if (fxxx[ii]-fxxx[ii-1])*(fxxx[ii+1]-fxxx[ii])<0:
                   xxleft.append(xxx[ii-1])
                   xxright.append(xxx[ii+1])
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
         #  print(xleft)
           xright=xxright[I]
         #  print(xright)
        else:
            y=xc
            break


    return y
#if (~isempty(xc))
#    y=xc;
#    break;
#else 
#    k=1;
#    for ii=2:1:length(xxx)-1
        
# %    secondd(ii)=(BB(ii)-BB(ii-1));
#        if((fxxx(ii)-fxxx(ii-1))*(fxxx(ii+1)-fxxx(ii))<0)
#           xxleft(k)=xxx(ii-1);
#           xxright(k)=xxx(ii+1);
#           m(k)=fxxx(ii);
#     %    break;
#           k=k+1;
#        end
 #   end
#    [~,I]=min(abs(m));
#    xleft=xxleft(I);
#    xright=xxright(I);
#end

#end

#end
