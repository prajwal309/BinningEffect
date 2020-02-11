import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from LineProfiles import *

WaveNumber = np.arange(950,1050,0.0025)

LineCenter = 1000
ScaleD = 10.0
ScaleL = 3.0
GamD = 1.0
GamL = 1.0

DopProf1 = PROFILE_DOPPLER(LineCenter, GamD, GamL, WaveNumber)
DopProf2 = PROFILE_DOPPLER(LineCenter, ScaleD*GamD, ScaleL*GamL, WaveNumber)

LorProf1 = PROFILE_LORENTZ(LineCenter, GamD, GamL, WaveNumber)
LorProf2 = PROFILE_DOPPLER(LineCenter, ScaleD*GamD), ScaleL*GamL, WaveNumber)

VoigtProf1 = PROFILE_VOIGT(LineCenter, GamD, GamL, WaveNumber)#+PROFILE_VOIGT(1050, 3.0, 3.0, WaveNumberVoigt)
VoigtProf2 = PROFILE_VOIGT(LineCenter, ScaleD*GamD, ScaleL*GamL, WaveNumber)#+PROFILE_VOIGT(1050, 1.0, 1.0, WaveNumberVoigt)




#Create line individually
ConvProf1 = (WaveNumber[1]-WaveNumber[0])*np.convolve(DopProf1, LorProf1, mode='same')
ConvProf2 = (WaveNumber[1]-WaveNumber[0])*np.convolve(DopProf2, LorProf2, mode='same')

TrialConv = (WaveNumber[1]-WaveNumber[0])*np.convolve(VoigtProf1,np.sqrt(DopProf2), mode='same')
#TrialConv = (WaveNumber[1]-WaveNumber[0])*np.convolve(TrialConv,LorProf2, mode='same')

plt.figure()
plt.plot(WaveNumber, TrialConv, color="navy", linestyle=":", label="Obtained Profile")
plt.plot(WaveNumber, VoigtProf2, "bo", label="Expected Profile")
plt.show()

#Taking the fast fourier transform
FFT_Dop1 = 1./len(WaveNumber)*np.fft.fft(DopProf1)
FFT_Dop2 = 1./len(WaveNumber)*np.fft.fft(DopProf2)


FFT_Lor1 = 1./len(WaveNumber)*np.fft.fft(LorProf1)
FFT_Lor2 = 1./len(WaveNumber)*np.fft.fft(LorProf2)

MultipliedFFT1 = np.fft.fftshift(np.fft.ifft(FFT_Dop1*FFT_Lor1))

plt.figure()
plt.plot(WaveNumber, VoigtProf1/max(VoigtProf1), "ko")
plt.plot(WaveNumber, ConvProf1/max(ConvProf1), "r-")
plt.plot(WaveNumber, MultipliedFFT1/max(MultipliedFFT1), "g+")
plt.show()

#RescaleFactor = (WaveNumber[1]-WaveNumber[0])
#VoigtConvolved1 = RescaleFactor*np.convolve(LD_Profile, LL_Profile1, mode='same')
#VoigtConvolved2 = RescaleFactor*np.convolve(LD_Profile, LL_Profile2, mode='same')
