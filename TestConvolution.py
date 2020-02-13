import numpy as np
import matplotlib.pyplot as plt
import emcee
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
    #A simple polynomial
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



def log_likelihood(theta, WaveNumber, StandardProfile, Sg0, Error):
    global LeastError
    GamD = 1.0#np.random.uniform(0,1)           #Doppler FWHM
    GamL = 1.0#np.random.uniform(0,1)           #Lorenztian FWHM

    FunctValue =Func2Conv3(theta, GamD, GamL, WaveNumber, Sg0)

    ActualProfile = PROFILE_VOIGT(WaveNumber, GamD, GamD, WaveNumber0)


    ConvolvedModel = (WaveNumber[1]-WaveNumber[0])*np.convolve(StandardProfile, FunctValue, mode="same")

    Difference = (ActualProfile - ConvolvedModel)
    TotalResidual = np.sum(Difference*Difference)
    ChiSq = -TotalResidual/1e-4
    if LeastError>TotalResidual:
        #plt.figure()
        #plt.plot(WaveNumber, ActualProfile, "ko")
        #plt.plot(WaveNumber, ConvolvedModel, "ko")
        #plt.show()

        LeastError = TotalResidual
        print(ChiSq)
        np.savetxt("BestParameters/BestParam.txt", theta)

    return ChiSq




PolyOrder = 5
StepSize = 0.01
Error = 0.005

#Generate Data
np.random.seed(10)

WaveNumber = np.arange(980,1020,StepSize)
WaveNumber0 = 1000
Profile = PROFILE_VOIGT(WaveNumber, 1.0, 1.0, WaveNumber0)

#plt.figure()
#plt.plot(WaveNumber, Profile, "ko")
#plt.show()


nWalkers = 4*PolyOrder
StartingGuessSelected = np.random.normal(0,0.01,(nWalkers,PolyOrder))

global LeastError
LeastError = np.inf


NSteps = 15000

#Initialize the transit model with batman
sampler = emcee.EnsembleSampler(nWalkers, PolyOrder, log_likelihood, args=[WaveNumber,Profile,WaveNumber0, Error], threads=1)
pos, prob, state = sampler.run_mcmc(StartingGuessSelected, NSteps)

#Start the en
