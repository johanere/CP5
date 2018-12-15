import numpy as np
from mpl_toolkits.mplot3d import Axes3D # This import has side effects required for the kwarg projection='3d' in the call to fig.add_subplot
import matplotlib.pyplot as plt
import random

A=np.loadtxt("exact")
n=np.size(A,0)
print n
 
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = y = np.arange(0, 1, 1.0/n)
X, Y = np.meshgrid(x, y)
zs = A
Z = zs.reshape(X.shape)
 
ax.plot_surface(X, Y, Z)
 
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u(x,t)')
 
plt.savefig("../../../Results/2dexact.pdf", bbox_inches='tight')

A=np.loadtxt("numerical")
n=np.size(A,0)
print n
 
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = y = np.arange(0, 1, 1.0/n)
X, Y = np.meshgrid(x, y)
zs = A
Z = zs.reshape(X.shape)
 
ax.plot_surface(X, Y, Z)
 
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u(x,t)')
 
plt.savefig("../../../Results/2dnumerical.pdf", bbox_inches='tight')

A=np.loadtxt("relerr")
n=np.size(A,0)
print n
 
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = y = np.arange(0, 1, 1.0/n)
X, Y = np.meshgrid(x, y)
zs = A
Z = zs.reshape(X.shape)
 
ax.plot_surface(X, Y, Z)
 
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u(x,t)')
 
plt.savefig("../../../Results/2drelerr.pdf", bbox_inches='tight')

