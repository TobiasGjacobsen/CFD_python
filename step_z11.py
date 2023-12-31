import numpy as np
import matplotlib.pyplot as plt
nx = 50; ny = 50; nt = 900; nit = 50;
dx = 2 / (nx - 1); dy = 2 / (ny - 1); x = np.linspace(0, 2, nx); y = np.linspace(0, 2, ny)
rho = 1; nu = .1; dt = .001
u = np.zeros((ny, nx)); v = np.zeros((ny, nx)); p = np.zeros((ny, nx)); b = np.zeros((ny, nx))

for n in range(nt):
    un = u.copy(); vn = v.copy()
    b[1:-1, 1:-1] =(rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) /
    (2 * dy)) - ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx)) ** 2 - 2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) /
    (2 * dy) * (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) - ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) ** 2))

    print(b)
    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy ** 2 + (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx ** 2) /
                             (2 * (dx ** 2 + dy ** 2)) - dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) * b[1:-1, 1:-1])

        p[:, -1] = p[:, -2]; p[0, :] = p[1, :]; p[:, 0] = p[:, 1]; p[-1, :] = 0
        # dp/dx = 0 at x = 2 , dp/dy = 0 at y = 0, dp/dx = 0 at x = 0, # p = 0 at y = 2

    u[1:-1, 1:-1] = (un[1:-1, 1:-1]-un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy *(un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                     dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +nu * (dt / dx ** 2 *
                     (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) + dt / dy ** 2 *
                     (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))


    u[0, :] = 0; u[:, 0] = 0; u[:, -1] = 0; u[-1, :] = 1  # set velocity on cavity lid equal to 1

    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx *(vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy *(vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -dt / (2 * rho * dy) * (p[2:, 1:-1] -
                     p[0:-2, 1:-1]) +nu * (dt / dx ** 2 *(vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +dt / dy ** 2 *
                     (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

    v[0, :] = 0; v[-1, :] = 0; v[:, 0] = 0; v[:, -1] = 0


plt.figure(figsize=(11, 7))
X, Y = np.meshgrid(x, y)
plt.contourf(X, Y, p, cmap='jet')
colorbar = plt.colorbar()
colorbar.set_label('Pressure [Pa]')
plt.contour(X, Y, p, cmap='jet')
plt.streamplot(X, Y, u, v)
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

ax = plt.axes(projection='3d')
X, Y = np.meshgrid(x, y)
surf = ax.plot_surface(X, Y, b[:], cmap='jet')
plt.show()

