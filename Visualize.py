import numpy as np
import matplotlib.pyplot as plt
from LineProfiles import PROFILE_VOIGT

def Func2Conv1(theta, GamD, GamL, WaveNumber, WaveNumber0):

    alpha = GamL/GamD
    Dist = np.abs(WaveNumber - WaveNumber0)


    #A simple polynomial
    FunctionVal = np.polyval(theta, alpha)*Dist

    ##
    return FunctionVal

def Func2Conv2(theta, GamD, GamL, WaveNumber, WaveNumber0):

    alpha = GamL/GamD
    Dist = (WaveNumber - WaveNumber0)*(WaveNumber - WaveNumber0)

    FunctionVal = np.exp(-np.polyval(theta, alpha)*Dist)

    ##
    return FunctionVal



def Func2Conv3(theta, GamD, GamL, WaveNumber, WaveNumber0):

    alpha = GamL/GamD
    Dist = np.abs((WaveNumber - WaveNumber0))
    #A simple polynomial
    FunctionVal = np.exp(-np.polyval(theta, GamL)*Dist)+np.exp(-np.polyval(theta, GamD)*Dist)
    ##
    return FunctionVal





theta = np.loadtxt("BestParameters/BestParam.txt")
print("The parameters are::", theta)

GamD = 1.0
GamL = 1.0
StepSize = 0.01
WaveNumber = np.arange(980,1020,StepSize)
WaveNumber0 = 1000

FunctValue =Func2Conv3(theta, GamD, GamL, WaveNumber, WaveNumber0)
ActualProfile = PROFILE_VOIGT(WaveNumber, GamD, GamD, WaveNumber0)
StandardProfile = PROFILE_VOIGT(WaveNumber, 1.0, 1.0, WaveNumber0)

ConvolvedModel = (WaveNumber[1]-WaveNumber[0])*np.convolve(StandardProfile, FunctValue, mode="same")


plt.figure(figsize=(12,10))
plt.subplot(211)
plt.plot(WaveNumber, FunctValue, "g-")
plt.subplot(212)
plt.plot(WaveNumber, ActualProfile, "ko")
plt.plot(WaveNumber, ConvolvedModel, "r-")
plt.show()
