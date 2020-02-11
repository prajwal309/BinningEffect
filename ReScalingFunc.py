import numpy as np
import matplotlib.pyplot as plt

def Gaussian(Freq, Freq0, Gamma_D, A):
 return 1./(np.sqrt(np.pi)*Gamma_D)*np.exp(-(Freq - Freq0)**2.0/Gamma_D**2.0)

def Lorentzian(Freq, Freq0, Gamma_L, A):
 return Gamma_L/((Freq-Freq0)**2.0+Gamma_L**2.0)



Freq = np.linspace(95,105,1000)



Nu_0 = 100
Gamma_D1 = 1.0
Gamma_D2 = 0.5

Gamma_L1 = 1.0
Gamma_L2 = 0.5


G1 = Gaussian(Freq,100, Gamma_D1, 1.0 )
G2 = Gaussian(Freq,100, Gamma_D2, 1.0)

L1 = Lorentzian(Freq,100, Gamma_L1, 1.0 )
L2 = Lorentzian(Freq,100, Gamma_L2, 1.0 )

H1 = np.convolve(G1, L2, mode='same')
H2 = np.convolve(G2, L1, mode='same')



plt.figure()
plt.subplot(211)
plt.plot(Freq, G1, "ko")
plt.plot(Freq, G2, "ro")
plt.subplot(212)
plt.plot(Freq, L1, "ko")
plt.plot(Freq, L2, "ro")
plt.show()


plt.figure(figsize=(12,7))
plt.subplot(211)
plt.plot(Freq, H1, "k-", label="Case 1")
plt.plot(Freq, H2+2.0, "r-", label="Case 2")
plt.legend()
plt.subplot(212)
plt.plot(Freq, H1-H2, "k-", label="Case 1")
plt.xlabel("Frequency")
plt.ylabel("Residual")
plt.show()



Area1 = np.sum(Gaussian1)
