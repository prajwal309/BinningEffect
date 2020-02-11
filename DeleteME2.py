import numpy as np
import matplotlib.pyplot as plt

x = np.random.uniform(0,1,10000)

for i in range(10):
    plt.figure()
    plt.plot(x, "ko-")
    plt.show()
    x = np.convolve(x,x,mode="same")
