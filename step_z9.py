import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

nx = 31; ny = 31; c = 1; dx = 2 / (nx - 1); dy = 2 / (ny - 1)

##initial conditions
p = np.zeros((ny, nx))  # create a XxY vector of 0's
##plotting aids
x = np.linspace(0, 2, nx); y = np.linspace(0, 1, ny)
##boundary conditions
p[:, 0] = 0  # p = 0 @ x = 0
p[:, -1] = y  # p = y @ x = 2
p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

#pd.set_option('display.max_rows', None)
#pd.set_option('display.max_columns', None)
data = p
df = pd.DataFrame(data)
print(" ")
print("initial pressure matrix:")
print("------------------------")
print(df)



plt.figure(figsize=(11, 7))
ax = plt.axes(projection='3d')
X, Y = np.meshgrid(x, y)
ax.plot_surface(X, Y, p[:], cmap='jet')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
surf = ax.plot_surface(X, Y, p[:], cmap='jet')
plt.colorbar(surf)
#plt.show()

l1norm = 1; l1norm_target=1e-4


while l1norm > 1e-4:
    pn = p.copy()
    p[1:-1, 1:-1] = ((dy ** 2 * (pn[1:-1, 2:] + pn[1:-1, 0:-2]) +
                      dx ** 2 * (pn[2:, 1:-1] + pn[0:-2, 1:-1])) /
                     (2 * (dx ** 2 + dy ** 2)))

    p[:, 0] = 0  # p = 0 @ x = 0
    p[:, -1] = y  # p = y @ x = 2
    p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
    p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1
    l1norm = np.sum(np.abs(p[:]) - np.abs(pn[:])) / np.sum(np.abs(pn[:]))

plt.figure(figsize=(11, 7), dpi=100)
ax = plt.axes(projection='3d')
X, Y = np.meshgrid(x, y)
ax.plot_surface(X, Y, p[:], cmap='jet')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
surf = ax.plot_surface(X, Y, p[:], cmap='jet')
plt.colorbar(surf)


data = p
df = pd.DataFrame(data)
print(" ")
print("Final pressure matrix:")
print("----------------------")
print(df)
plt.show()