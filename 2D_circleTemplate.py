import turtle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def circle(radius, centerx, centery, nx):
    theta = np.linspace(0, 2*np.pi, nx)
    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)
    return x, y

nx = 15
ny = 15
radius = 0.5
L = 2
H = 2
x = np.linspace(0, L, nx)
y = np.linspace(0, H, ny)
dy = H / (ny - 1)
dx = L / (nx - 1)

center = (1, 1)
circ=circle(radius, center[0], center[1], nx)

u = np.ones((ny,ny))
data = pd.DataFrame(u)
#print(data)
#plt.plot(circ[0], circ[1])
#plt.show()
gridI=np.ones(nx)
gridJ=np.ones(ny)
xc = circ[0]
yc = circ[1]
for i in range(nx):
    gridI[i] = xc[i]/dx
    gridJ[i] = yc[i]/dy
    print(gridJ, gridI)
    print(x)
    u[int(gridJ[i]), int(gridI[i])]=0

plt.plot(gridI, gridJ)
plt.show()

data = pd.DataFrame(u)
print(data)

#gridJ * dx = Ly -> Ly/dy =gridJ
#gridI * dy = Lx -> Lx/dx =GridI

'''
working boundary conditions for the wall of the cylinder 

    for it in range(nx):
        gridI[it] = xc[it] / dx
        gridJ[it] = yc[it] / dy
        u[int(gridJ[it]), int(gridI[it])] = 0
        v[int(gridJ[it]), int(gridI[it])] = 0

    for vert in np.linspace(0,2*np.pi, nx):
        xglob = center[0] + radius*np.cos(vert)
        yglob = center[1] + radius*np.sin(vert)
        u[int(center[1]/dy):int(yglob/dy),int(xglob/dx)]=0
        v[int(center[1] / dy):int(yglob / dy), int(xglob / dx)] = 0
        xglob = center[0] - radius*np.cos(vert)
        yglob = center[1] - radius*np.sin(vert)
        u[int(yglob/dy):int(center[1]/dy),int(xglob/dx)]=0
        v[int(yglob / dy):int(center[1] / dy), int(xglob / dx)] = 0

'''

