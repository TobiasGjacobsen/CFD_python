import numpy as np
import matplotlib.pyplot as plt
grid_inc = 1

simtime = 1
nx = 50*grid_inc; ny = 50*grid_inc; dt = .001; nt =int(simtime/dt); nit = 50; c = 1
dx = 2 / (nx - 1); dy = 2 / (ny - 1); x = np.linspace(0, 2, nx); y = np.linspace(0, 2, ny)
rho = 1; nu = .1;
u = np.zeros((ny, nx)); v = np.zeros((ny, nx)); p = np.zeros((ny, nx)); b = np.zeros((ny, nx))

if grid_inc > 1:
    dt = .001 / (grid_inc * 2)
    print("modifying timestep")
    print("dt: ", dt)
    nt = int(simtime / dt)
print("physical time: ", dt*nt)

def build_up_b(b, rho, dt, u, v, dx, dy):
    b[1:-1, 1:-1] = (rho * (1 / dt *((u[1:-1, 2:] - u[1:-1, 0:-2]) /(2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) /
    (2 * dy)) - ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx)) ** 2 -2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) /
    (2 * dy) * (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))-((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) ** 2))
    return b

def pressure_poisson(p, dx, dy, b):
    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy ** 2 + (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx ** 2)/
        (2 * (dx ** 2 + dy ** 2)) - dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) * b[1:-1, 1:-1])

        p[:, -1] = p[:, -2]  # dp/dx = 0 at x = 2
        p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
        p[:, 0] = p[:, 1]  # dp/dx = 0 at x = 0
        p[-1, :] = 0  # p = 0 at y = 2

    return p

def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu):
    b = np.zeros((ny, nx))
    for n in range(nt):
        un = u.copy(); vn = v.copy()
        b = build_up_b(b, rho, dt, u, v, dx, dy)
        p = pressure_poisson(p, dx, dy, b)

        u[1:-1, 1:-1] = (un[1:-1, 1:-1]-un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
        vn[1:-1, 1:-1] * dt / dy *(un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
        dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +nu * (dt / dx ** 2 *
        (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) + dt / dy ** 2 *
        (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

        v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx *(vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
        vn[1:-1, 1:-1] * dt / dy *(vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
        dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +nu * (dt / dx ** 2 *
        (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +dt / dy ** 2 *
        (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

        u[0, :] = 0; u[:, 0] = 0; u[:, -1] = 0; u[-1, :] = 1  # set velocity on cavity lid equal to 1
        v[0, :] = 0; v[-1, :] = 0; v[:, 0] = 0; v[:, -1] = 0
        print("running")

    return u, v, p
u, v, p = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu)

plt.figure(figsize=(12,10))
X, Y = np.meshgrid(x, y)
plt.contourf(X, Y, p, cmap='jet')
colorbar = plt.colorbar()
colorbar.set_label('Pressure [Pa]')
plt.contour(X, Y, p, cmap='jet')
plt.streamplot(X, Y, u, v)
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

plt.figure(figsize=(12, 10))
X, Y = np.meshgrid(x, y)
plt.streamplot(X, Y, u, v, linewidth=1.5, color=np.sqrt(u**2+v**2), cmap=plt.cm.jet, density=3.0)
colorbar = plt.colorbar()
colorbar.set_label('velocity')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()