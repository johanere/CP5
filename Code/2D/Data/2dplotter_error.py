import numpy as np
from mpl_toolkits.mplot3d import Axes3D # This import has side effects required for the kwarg projection='3d' in the call to fig.add_subplot
import matplotlib.pyplot as plt
import random

def tex(x):
    exp = int(np.floor(np.log10(abs(x))))
    number=x/10**exp
    return   "$%.2g \cdot 10^{%s}$"%(number,exp)

outdata="Avr_relerr"

f=open(outdata)
lines = f.readlines()
N=len(lines)

relerr=np.zeros(N)
dtfactor=np.zeros(N)

for i in range(0,N):
    hold=lines[i].split()
    relerr[N-1-i]=float(hold[0])
    dtfactor[N-1-i]=float(hold[1])
  


print relerr
print dtfactor

  
plt.plot(dtfactor,relerr)
plt.show()