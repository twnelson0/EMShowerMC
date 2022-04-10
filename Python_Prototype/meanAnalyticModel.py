import numpy as np
import matplotlib.pyplot as plt

def E(E0,t):
    return E0/2**t

def dEdt(E0,t):
    return E0/2**t *np.log(2)

tArr = np.linspace(0,25)
E0 = 50000

plt.plot(tArr,1/E0*dEdt(E0,tArr))
plt.show()
