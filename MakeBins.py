import numpy as np
import matplotlib.pyplot as plt
from LineProfiles import *

#Assume a molecular mass

#Calculate the gamma D


Gam0 = 0.
GamD = 0.1

GamV = 0.5346*Gam0 + np.sqrt(0.2166*GamD**2)

LineCenter = 1000.0

#How far is one going
Extent = 1500.0*GamV

StepLevels = np.arange(-5,5.,0.01)
AllStepSizes= GamV/2.0**StepLevels

AreaValues = np.zeros(len(AllStepSizes))

Tolerance = 1e-5
PrintFlag = True

for Counter, StepValue in enumerate(AllStepSizes):
    WaveNumber = np.arange(LineCenter-Extent,LineCenter+Extent,StepValue)

    #Generate theWaveNumber profiling
    Profile = PROFILE_VOIGT(LineCenter, GamD, Gam0, WaveNumber)
    Area = np.trapz(Profile, WaveNumber)
    if PrintFlag and Counter>5 and max(abs(AreaValues[Counter-5:Counter]-1.0))<Tolerance:
        StoreValue = StepValue
        print(StepValue/GamV)
        PrintFlag = False
    AreaValues[Counter]=Area

plt.figure(figsize=(12,8))
plt.plot(AllStepSizes/GamV, AreaValues, "ko")
plt.axvline(x=StoreValue/GamV, color="blue")
plt.axhline(y=1, color="red")
plt.xscale('log')
plt.yscale('log')
#plt.gca().invert_xaxis()
plt.ylabel("Area Under Curve")
plt.xlabel("StepSize/$\\gamma_V$")
plt.show()
