# -----------------------------------------------------------------------------
# ----------------------------------------------------------------------------- 
"""
Contributor:    Qian Zhang, Northwestern University
                qian.zhang_jennifer@northwestern.edu
Created:   Apr 2020 
Modified:  ... ...

Description
-----------
Build a GUI wrapper to explore the distribution of oxygen parial pressure in two-layer electrolyte.
"""
try:
    # This will work in Python 2.7
    import Tkinter
except ImportError:
    # This will work in Python 3.5
    import tkinter as Tkinter

# -----------------------------------------------------------------------------
# To use matplotlib, the author must use the TkAgg backend, or none of this will
# work and a long string of inexplicable error messages will ensue.
# ----------------------------------------------------------------------------- 
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
from PFHYBRID import PFHYBRID
from OxygenParT import OxygenParT
from PFHYBRIDinitial import PFHYBRIDinitial

# Define a bold font:
BOLD = ('Courier', '15', 'bold')

# Create main application window.
root = Tkinter.Tk()

# Create a text box explaining the application.
greeting = Tkinter.Label(text="PO2 distribution in the electrolyte", font=BOLD)
greeting.pack(side='top')

# Create a frame for variable names and entry boxes for their values.
frame = Tkinter.Frame(root)
frame.pack(side='top')

# Variables for the calculation, and default values.
Temp = Tkinter.StringVar()
Temp.set('800.0')

CuDe = Tkinter.StringVar()
CuDe.set('-0.8')

TYSZ = Tkinter.StringVar()
TYSZ.set('10.0')

TGDC = Tkinter.StringVar()
TGDC.set('2.5')

ComHydr = Tkinter.StringVar()
ComHydr.set('50.0')

RpH = Tkinter.StringVar()
RpH.set('0.1')

RpO = Tkinter.StringVar()
RpO.set('0.1')

Symindex = Tkinter.StringVar()
Symindex.set('1.0')

# Create text boxes and entry boxes for the variables.
# Use grid geometry manager instead of packing the entries in.
row_counter = 0
aa_text = Tkinter.Label(frame, text='Temperature (C):') 
aa_text.grid(row=row_counter, column=0)

aa_entry = Tkinter.Entry(frame, width=8, textvariable=Temp)
aa_entry.grid(row=row_counter, column=1)

row_counter += 1
fa_text = Tkinter.Label(frame, text='Current Density (A/cm^2):') 
fa_text.grid(row=row_counter, column=0)

fa_entry = Tkinter.Entry(frame, width=8, textvariable=CuDe)
fa_entry.grid(row=row_counter, column=1)

row_counter += 1
ab_text = Tkinter.Label(frame, text='Thickness YSZ (Micron):') 
ab_text.grid(row=row_counter, column=0)

ab_entry = Tkinter.Entry(frame, width=8, textvariable=TYSZ)
ab_entry.grid(row=row_counter, column=1)

row_counter += 1
fb_text = Tkinter.Label(frame, text='Thickness GDC (Micron):') 
fb_text.grid(row=row_counter, column=0)

fb_entry = Tkinter.Entry(frame, width=8, textvariable=TGDC)
fb_entry.grid(row=row_counter, column=1)

row_counter += 1
dp_text = Tkinter.Label(frame, text='Hydrogen Composition (between 0 and 100):') 
dp_text.grid(row=row_counter, column=0)

dp_entry = Tkinter.Entry(frame, width=8, textvariable=ComHydr)
dp_entry.grid(row=row_counter, column=1)

row_counter += 1
rph_text = Tkinter.Label(frame, text='Polarization Resistance at Hydrogen Side (Omega cm^2):') 
rph_text.grid(row=row_counter, column=0)

rph_entry = Tkinter.Entry(frame, width=8, textvariable=RpH)
rph_entry.grid(row=row_counter, column=1)

row_counter += 1
rpo_text = Tkinter.Label(frame, text='Polarization Resistance at Oxygen Side (Omega cm^2):') 
rpo_text.grid(row=row_counter, column=0)

rpo_entry = Tkinter.Entry(frame, width=8, textvariable=RpO)
rpo_entry.grid(row=row_counter, column=1)

row_counter += 1
rpo_text = Tkinter.Label(frame, text='Symmetric cell? (YSZ electrolyte):') 
rpo_text.grid(row=row_counter, column=0)

rpo_entry = Tkinter.Entry(frame, width=8, textvariable=Symindex)
rpo_entry.grid(row=row_counter, column=1)

# Define a function to create the desired plot.
def make_plot(event=None):
    # Get these variables from outside the function, and update them.
    global Temp, CuDe, TYSZ, TGDC, ComHydr, RpH, RpO, Symindex

    # Convert StringVar data to numerical data.
    aa = float(Temp.get())
    fa = float(CuDe.get())
    ab = float(TYSZ.get())
    fb = float(TGDC.get())
    phi = float(ComHydr.get())
    Rph = float(RpH.get())
    Rpo = float(RpO.get())
    Symindexf = float(Symindex.get())
    
   
    # Define the range of the plot.
    #t_min = -10
    #t_max = 10
    #dt = 0.01
    #t = np.arange(t_min, t_max+dt, dt)

    # Create the two waves and find the combined intensity.
    #waveA = aa * np.cos(fa * t)
    #waveB = ab * np.cos(fb * t + phi)
    #intensity = (waveA + waveB)**2
    L_length=(ab+fb)*1e-4
    PO2EH=OxygenParT(aa,phi,100-phi)
    x_portion=fb/(ab+fb)
    #err=0.0

    #if (err==0.0):
    #   master = Tkinter.Tk()
    #   whatever_you_do = "Please wait ... ... ."
    #   msg = Tkinter.Message(master, text = whatever_you_do)
    #   msg.config(bg='lightgreen', font=('times', 24, 'italic'))
    #   msg.pack()
    #
   # if x_portion<=0.5:
    if (Symindexf == 0.0):
       if (aa<=950 and x_portion<=0.5):
           methodRK="RK23"
           xp, po2xp, err=PFHYBRID(aa,fa,L_length,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,methodRK)
    #elif aa<=950:
    #   xp, po2xp, err=PFHYBRID(aa,fa,L_length,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,"RK23")
       elif (aa>950 and x_portion<=0.5):
            methodRK="RK45"
            xp, po2xp, err=PFHYBRID(aa,fa,L_length,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,methodRK)
       elif (aa<=950 and x_portion>0.5):
            methodRK="RK23"
          #  err, initial=PFHYBRIDinitial(aa,fa,L_length,0.1*Rph,0.1*Rpo,PO2EH,0.2,0.5,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,"RK23")
            xp, po2xp, err=PFHYBRID(aa,fa,L_length,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,methodRK)
       elif (aa>950 and x_portion>0.5):
            methodRK="RK45"
           # err, initial=PFHYBRIDinitial(aa,fa,L_length,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,"RK45")
            xp, po2xp, err=PFHYBRID(aa,fa,L_length,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,0,1,1,1,1,1, 1e-3, 1e-6, 1e-3,methodRK)      
      #elif x_portion>0.5:
          # if aa<=950:
          #    err, initial=PFHYBRIDinitial(aa,fa,L_length,Rph,Rpo,PO2EH,0.2,0.5,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,"RK23")
          #    xp, po2xp, err=PFHYBRID(aa,fa,L_length,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,initial,1,1,0,0,1, 1e-3, 1e-6, 1e-3,"RK23")
          #    methodRK="RK23"
    #elif aa<=950:
    #   xp, po2xp, err=PFHYBRID(aa,fa,L_length,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,"RK23")
           # elif aa>950:
           #      err, initial=PFHYBRIDinitial(aa,fa,L_length,Rph,Rpo,PO2EH,0.2,0.5,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,"RK45")
           #      xp, po2xp, err=PFHYBRID(aa,fa,L_length,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,initial,1,1,0,0,1, 1e-3, 1e-6, 1e-3,"RK45")
           #      methodRK="RK45"
    else:
        err=0.0

    faformer=fa
    L_lengthformer=L_length
        
    while err==1:
        master = Tkinter.Tk()
        whatever_you_do = "It may take longer time to do the calculation."
        msg = Tkinter.Message(master, text = whatever_you_do)
        msg.config(bg='lightgreen', font=('times', 24, 'italic'))
        msg.pack() 
        fa=fa-np.sign(fa)*0.2
        err, initial=PFHYBRIDinitial(aa,fa,L_lengthformer,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,methodRK)
        if err==0:
          xp, po2xp, err=PFHYBRID(aa,faformer,L_lengthformer,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,initial,1,1,0,0,1, 1e-3, 1e-6, 1e-3,methodRK)  
          if err==0:
              break
          else:
              continue
        else:
          L_length=L_length*0.75
          err, initial=PFHYBRIDinitial(aa,faformer,L_length,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,methodRK)
          if err==0:
            xp, po2xp, err=PFHYBRID(aa,faformer,L_lengthformer,Rph,Rpo,PO2EH,0.2,x_portion,2e+4,initial,1,1,0,0,1, 1e-3, 1e-6, 1e-3,methodRK)  
            if err==0:
                break 
            else:
                continue
              
        #else:

    if (Symindexf == 1.0 and err==0.0):
        methodRK="RK23"
          #  err, initial=PFHYBRIDinitial(aa,fa,L_length,0.1*Rph,0.1*Rpo,PO2EH,0.2,0.5,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,"RK23")
        xp, po2xp, err=PFHYBRID(aa,fa,L_length,Rph,Rpo,PO2EH,PO2EH,-1,2e+4,0,1,1,1,1,1, 1e-3, 1e-6, 1e-3,methodRK)
    elif (Symindexf == 2.0 and err==0.0):
        methodRK="RK23"
          #  err, initial=PFHYBRIDinitial(aa,fa,L_length,0.1*Rph,0.1*Rpo,PO2EH,0.2,0.5,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,"RK23")
        xp, po2xp, err=PFHYBRID(aa,fa,L_length,Rph,Rpo,0.2,0.2,-1,2e+4,0,1,1,1,0,1, 1e-3, 1e-6, 1e-3,methodRK)  
    else:
        master = Tkinter.Tk()
        whatever_you_do = "Please double check your input parameters and make sure your cell is a normal one."
        msg = Tkinter.Message(master, text = whatever_you_do)
        msg.config(bg='lightgreen', font=('times', 24, 'italic'))
        msg.pack()       
    #print('err=',err)
    # Create the plot.
    plt.figure()
    plt.plot(xp, po2xp, lw=3)
    plt.title('PO2 distribution')
    plt.xlabel('x/L')
    plt.ylabel('PO2 (atm)')
    plt.yscale("log")
    plt.show()


# Add a button to create the plot.
MakePlot = Tkinter.Button(root, command=make_plot, text="Create PO2 Plot")
MakePlot.pack(side='bottom', fill='both')


# Allow pressing <Return> to create plot.
root.bind('<Return>', make_plot)

# Allow pressing <Esc> to close the window.
root.bind('<Escape>', root.destroy)

# Activate the window.
root.mainloop()

