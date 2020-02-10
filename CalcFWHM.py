import numpy as np
import matplotlib.pyplot as plt

#Define the constants
cBolts = 1.380648813E-16 # erg/K, CGS
cc = 2.99792458e10 # cm/s, CGS
hh = 6.626196e-27 # erg*s, CGS


MoleculeName = ["N2", "H2O", "CH4", "CO2", "CO", "H2"]
Mass = [28.006148, 18.010565, 16.0313, 43.98983, 27.994915, 2.01565]


Temp = 100  #The temperature is 100 K
LineCenter = 2000   #per cm ---> convert from microns
                    #333.333 for 30 microns
for Name, m in zip(MoleculeName, Mass):
    GammaD = np.sqrt(2*cBolts*Temp*np.log(2)/m/cc**2)*LineCenter
    print(Name, GammaD)
