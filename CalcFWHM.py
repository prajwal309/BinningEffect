import numpy as np
import matplotlib.pyplot as plt


def GetWaveNumbers(StartWavelength=300, EndWavelength=30000, Resolution=100000):
        '''
        Returns the wavelengths corresponding to the resolution and in units
        of cm.
        '''
        WaveLengthValues = []

        #Converting values to
        StartWavelength = StartWavelength*1.0e-7       #nm to cm
        EndWavelength = EndWavelength*1.0e-7           #nm to cm
        WaveLengthValues = [StartWavelength]

        while WaveLengthValues[-1]<EndWavelength:
            WaveLengthValues.append(WaveLengthValues[-1]+WaveLengthValues[-1]/Resolution)
        WaveLengthValues = np.array(WaveLengthValues)
        WaveNumberRange = 1./WaveLengthValues

        WaveNumberRange = np.array(sorted(WaveNumberRange))
        return WaveLengthValues, WaveNumberRange

Wavelength, wavenumber = GetWaveNumbers(Resolution = 1000)
diff = np.diff(Wavelength)
Resolution = Wavelength[:-1]/diff



#Define the constants
cBolts = 1.380648813E-16 # erg/K, CGS
cc = 2.99792458e10 # cm/s, CGS
hh = 6.626196e-27 # erg*s, CGS
cMassMol = 1.66053873e-27 # Molecular mass equivalent in kilograms

MoleculeName = ["N2", "H2O", "CH4", "CO2", "CO", "H2"]
Mass = np.array([28.006148, 18.010565, 16.0313, 43.98983, 27.994915, 2.01565])*cMassMol*1000.0


Temp = 100  #Kelvins
WaveLength = np.arange(0.3, 30.0, 0.01)*1e-4
WaveNumber = 1./WaveLength

print("The minimum value of wavenumber is ::", min(WaveNumber))
print("The maximum value of wavenumber is ::", max(WaveNumber))

LineCenter = 333.33     #2000   #per cm ---> convert from microns
                        #333.333 for 30 microns

for Name, m in zip(MoleculeName, Mass):
    GammaD = np.sqrt(2*cBolts*Temp*np.log(2)/m/cc**2)*LineCenter
    Resolution = WaveLength/GammaD
    print(Resolution)
    input("Wait here...")
    plt.figure()
    plt.plot(WaveLength, Resolution, "k-")
    plt.xscale('log')
    plt.show()

    print(Name, GammaD)
