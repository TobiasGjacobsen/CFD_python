import matplotlib.pyplot as plt
import numpy as np

def circle(radius, centerx, centery, nx):
    theta = np.linspace(0, 2*np.pi, nx)
    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)
    return x, y

nx = 40
radius = 1.0
center = (10.0, 5.0)
circ=circle(radius, center[0], center[1], nx)

#gridJ * dx = Ly -> Ly/dy =gridJ
#gridI * dy = Lx -> Lx/dx =GridI

# Create the plot
plt.figure(figsize=(10, 10))
plt.plot(circ[0],circ[1])
plt.show()
