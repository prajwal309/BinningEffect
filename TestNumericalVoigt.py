import numpy as np
import matplotlib.pyplot as plt
from LineProfiles import *
import hapi
import time

def NUMERICAL_VOIGT(Sg, GamD, GamL, Sg0):
    d = (GamL-GamD)/(GamL+GamD)
    alpha = 0.18121
    beta = 0.023665*np.exp(0.6*d)+0.00418*np.exp(-1.9*d)
    cL = 0.6818817+0.6129331*d-0.1838439*d*d-0.1156844*d*d*d
    cG = 0.3246017-0.6182531*d+0.1768139*d*d+0.1210944*d*d*d
    GamV = (GamL + GamD)*(1- alpha*(1-d*d) - beta*np.sin(np.pi*d))

    Profile = cL*1./np.pi*GamV/((Sg - Sg0)*(Sg - Sg0)+GamV*GamV) + \
              cG*np.sqrt(np.log(2))/(np.sqrt(np.pi)*GamV)*np.exp(-np.log(2)*(Sg-Sg0)*(Sg-Sg0)/(GamV*GamV))
    return Profile

    

Sg = np.arange(990,1010,0.01)
GamD = 0.1
GamL = 1.0
LineCenter = 1000

NTrials = 100
StartTime = time.time()
for i in range(NTrials):
    HAPI_PROFILE = hapi.PROFILE_VOIGT(LineCenter, GamD, GamL, Sg)[0]
print("Time taken::", time.time() - StartTime)

StartTime = time.time()
for i in range(NTrials):
    VoigtProfile = PROFILE_VOIGT(LineCenter, GamD, GamL, Sg)
print("Time taken::", time.time() - StartTime)

StartTime = time.time()
for i in range(NTrials):
    NumProfile = NUMERICAL_VOIGT(LineCenter, GamD, GamL, Sg)
print("Time taken::", time.time() - StartTime)

plt.figure(figsize=(12,8))
plt.plot(Sg, VoigtProfile, "ko", label="My Current method")
plt.plot(Sg, NumProfile, "r+", label="New Numerical")
plt.plot(Sg, HAPI_PROFILE, "g-", label="HAPI Profile")
plt.xlabel("Wavenumber", fontsize=20)
plt.ylabel("Profile Strength", fontsize=20)
plt.xlim([996.5,1003.5])
plt.legend(loc=1)
plt.savefig("DifferentVoigtProfile.png")
plt.show()
