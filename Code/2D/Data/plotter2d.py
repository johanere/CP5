import numpy as np
from mpl_toolkits.mplot3d import Axes3D # This import has side effects required for the kwarg projection='3d' in the call to fig.add_subplot
import matplotlib.pyplot as plt
import random

def fun(x, y):
    t=1
    sum=0
    mu=x*np.pi
    nu=y*np.pi
    for n in range(1,400,2):
        for m in range(1,400,2):
            sum +=1/(m*n)*np.sin(m*x)*np.sin(n*y)*np.exp(-(m**2+n**2)*t)
        for m in range(2,400,2):
            sum -=1/(m*n)*np.sin(m*x)*np.sin(n*y)*np.exp(-(m**2+n**2)*t)
    for n in range(2,400,2):
        for m in range(1,400,2):
            sum -=1/(m*n)*np.sin(m*x)*np.sin(n*y)*np.exp(-(m**2+n**2)*t)
        for m in range(2,400,2):
            sum +=1/(m*n)*np.sin(m*x)*np.sin(n*y)*np.exp(-(m**2+n**2)*t)
    sum*4
    return sum+x+y

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = y = np.arange(0, np.pi, np.pi/10.)
X, Y = np.meshgrid(x, y)
zs = np.array([fun(x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
Z = zs.reshape(X.shape)

ax.plot_surface(X, Y, Z)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('temperature')

plt.show()