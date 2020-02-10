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
Gam0 = 0.1
GamD = 0.1

GamV = 0.5346*Gam0 + np.sqrt(0.2166*Gam0**2+GamD**2.0)


StepValue = 0.001

AreaValues = np.zeros(500)
LineCenter  = 1000.0

Factor = range(1,len(AreaValues)+1)


for Counter, Val in enumerate(Factor):
    #print("The step value is::", StepValue)
    Extent = Val*GamV
    WaveNumber = np.arange(LineCenter-Extent,LineCenter+Extent,StepValue)

    #Generate theWaveNumber profiling
    Profile = PROFILE_VOIGT(LineCenter, GamD, Gam0, WaveNumber)



    Area = np.trapz(Profile, WaveNumber)
    AreaValues[Counter]=Area

plt.figure(figsize=(12,8))
plt.axhline(y=1, color="red", linewidth=7, alpha=0.9)
plt.plot(Factor, AreaValues, "k-")
#plt.xscale('log')
plt.yscale('log')
#plt.gca().invert_xaxis()
plt.ylabel("Numerically Computed Area", fontsize=20 )
plt.xlabel("FWHM ($\\gamma_V$)", fontsize=20)
plt.tick_params(which="both", direction="in")
plt.tight_layout()
plt.savefig("FWHM.png")
plt.close('all')
plt.show()
