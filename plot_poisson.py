import numpy as np
from numpy import sin as sin
from numpy import cos as cos
from numpy import pi as pi
from numpy import exp as exp
#import poisson as ps
import math as mt
import scipy  
import matplotlib 
from matplotlib import pyplot as plt, cm
from mpl_toolkits.mplot3d import Axes3D



######## lid driven cavity ####### 
Lx = 50
Ly = 50
Nx = 50
Ny = 50 

dx = Lx / Nx
dy = Ly / Ny

x = np.linspace(0+dx/2.,Lx-dx/2.,Nx)
y = np.linspace(0+dy/2.,Ly-dy/2.,Ny)
################### 2n order coefficients for derivatives #############
p = np.zeros((Nx,   Ny))
p = np.loadtxt("p.dat")
print(p.shape)


##########################################
X, Y = np.meshgrid(x,y)

fig , (ax1)= plt.subplots(1,1) 
#ax1.figure(figsize=(11,7), dpi=100)
# plotting the pressure field as a contour
#ax1.contourf(X, Y, p_jac, alpha=0.5, cmap=cm.viridis)
#ax2.contourf(X, Y, p_gs, alpha=0.5, cmap=cm.viridis)
#ax1.contourf(X, Y, p,  levels=[0,2,4,6,8,10])

CS = ax1.contourf(X, Y, p, 50, cmap=plt.cm.bone) #, origin=origin)

# Note that in the following, we explicitly pass in a subset of
# the contour levels used for the filled contours.  Alternatively,
# We could pass in additional levels to provide extra resolution,
# or leave out the levels kwarg to use all of the original levels.

CS2 = plt.contour(CS, levels=CS.levels[::])#, colors='r')#, origin=origin)

plt.title('p')
plt.xlabel('x')
plt.ylabel('y')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel('verbosity coefficient')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)

# plotting the pressure field outlines
#ax1.set_xlabel('X')
#ax1.set_ylabel('Y')

plt.show()
