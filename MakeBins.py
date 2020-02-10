import numpy as np
import matplotlib.pyplot as plt
from LineProfiles import *

#Assume a molecular mass
import matplotlib as mpl
mpl.rc('font',**{'sans-serif':['Helvetica'], 'size':15,'weight':'bold'})
mpl.rc('axes',**{'labelweight':'bold', 'linewidth':1.5})
mpl.rc('ytick',**{'major.pad':22, 'color':'k'})
mpl.rc('xtick',**{'major.pad':10,})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
mpl.rc('text', **{'usetex':True})
mpl.rc('text.latex',preamble=r'\usepackage{cmbright},\usepackage{relsize},'+r'\usepackage{upgreek}, \usepackage{amsmath}')
mpl.rc('contour', **{'negative_linestyle':'solid'})


#Calculate the gamma D
Gam0 = 0.05
GamD = 0.1

GamV = 0.5346*Gam0 + np.sqrt(0.2166*Gam0**2+GamD**2.0)

LineCenter = 300.0

#How far is one going
Extent = 100.0*GamV

#StepLevels = np.arange(-8.0,3.,0.01)
StepLevels = np.arange(-4.0,1.5,0.0005)
AllStepSizes= GamV/2.0**StepLevels

AreaValues = np.zeros(len(AllStepSizes))

Tolerance = 0.005
PrintFlag = True



for Counter, StepValue in enumerate(AllStepSizes):
    #print("The step value is::", StepValue)
    WaveNumber = np.arange(LineCenter-Extent,LineCenter+Extent,StepValue)

    #Generate theWaveNumber profiling
    Profile = PROFILE_VOIGT(LineCenter, GamD, Gam0, WaveNumber)
    Area = np.trapz(Profile, WaveNumber)
    if PrintFlag and Counter>50 and max(abs(AreaValues[Counter-50:Counter]-1.0))<Tolerance:
        StoreValue = StepValue
        print(GamV/StoreValue)
        PrintFlag = False

        plt.figure()
        plt.plot(WaveNumber, Profile, "ko")
        plt.show()

    AreaValues[Counter]=Area

plt.figure(figsize=(12,8))
plt.axhline(y=1, color="red", linewidth=7, alpha=0.9)
plt.plot(GamV/AllStepSizes, AreaValues, "k-")
plt.axvline(x=GamV/StoreValue, color="blue")
plt.text(GamV/StoreValue*1.2,1.5,"Accuracy ~99.5\%",color="blue", fontsize=20)
#plt.xscale('log')
#plt.yscale('log')
#plt.gca().invert_xaxis()
plt.ylabel("Numerically Computed Area", fontsize=20 )
plt.xlabel("$\\gamma_V$/StepSize", fontsize=20)
plt.tick_params(which="both", direction="in")
plt.tight_layout()
plt.savefig("OptimalResolution.png")
plt.close('all')
plt.show()
