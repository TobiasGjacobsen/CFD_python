import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def build_up_b(rho, dt, dx, dy, u, v):
    b = np.zeros_like(u)
    b[1:-1, 1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) +
                                      (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                            ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx)) ** 2 -
                            2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                                 (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) -
                            ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) ** 2))

    # Periodic BC Pressure @ x = 2
    b[1:-1, -1] = (rho * (1 / dt * ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx) +
                                    (v[2:, -1] - v[0:-2, -1]) / (2 * dy)) -
                          ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx)) ** 2 -
                          2 * ((u[2:, -1] - u[0:-2, -1]) / (2 * dy) *
                               (v[1:-1, 0] - v[1:-1, -2]) / (2 * dx)) -
                          ((v[2:, -1] - v[0:-2, -1]) / (2 * dy)) ** 2))

    # Periodic BC Pressure @ x = 0
    b[1:-1, 0] = (rho * (1 / dt * ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx) +
                                   (v[2:, 0] - v[0:-2, 0]) / (2 * dy)) -
                         ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx)) ** 2 -
                         2 * ((u[2:, 0] - u[0:-2, 0]) / (2 * dy) *
                              (v[1:-1, 1] - v[1:-1, -1]) / (2 * dx)) -
                         ((v[2:, 0] - v[0:-2, 0]) / (2 * dy)) ** 2))

    return b


def pressure_poisson_periodic(p, dx, dy, elbow):
    pn = np.empty_like(p)

    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy ** 2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx ** 2) /
                         (2 * (dx ** 2 + dy ** 2)) -
                         dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) * b[1:-1, 1:-1])

        # Periodic BC Pressure @ x = 2
        p[1:-1, -1] = (((pn[1:-1, 0] + pn[1:-1, -2]) * dy ** 2 +
                        (pn[2:, -1] + pn[0:-2, -1]) * dx ** 2) /
                       (2 * (dx ** 2 + dy ** 2)) -
                       dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) * b[1:-1, -1])

        # Periodic BC Pressure @ x = 0
        p[1:-1, 0] = (((pn[1:-1, 1] + pn[1:-1, -1]) * dy ** 2 +
                       (pn[2:, 0] + pn[0:-2, 0]) * dx ** 2) /
                      (2 * (dx ** 2 + dy ** 2)) -
                      dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) * b[1:-1, 0])

        # Wall boundary conditions, pressure
        p[-1, :] = p[-2, :]  # dp/dy = 0 at y = 2
        p[0, :] = p[1, :]  # dp/dy = 0 at y = 0



    return p

##variable declarations
L = 6
H = 2
nx = 41*3
ny = 41*3
nt = 250
nit = 50
c = 1
dx = L / (nx - 1)
dy = H / (ny - 1)
x = np.linspace(0, L, nx)
y = np.linspace(0, H, ny)
X, Y = np.meshgrid(x, y)


##physical variables
rho = 1
nu = .1
F = 10
dt = .001
elbow = 50
albue = 40

#initial conditions
u = np.zeros((ny, nx))
un = np.zeros((ny, nx))

v = np.zeros((ny, nx))
vn = np.zeros((ny, nx))

p = np.ones((ny, nx))
pn = np.ones((ny, nx))

b = np.zeros((ny, nx))

udiff = 1
stepcount = 0

while udiff > .001:
    un = u.copy()
    vn = v.copy()

    b = build_up_b(rho, dt, dx, dy, u, v)
    p = pressure_poisson_periodic(p, dx, dy, elbow)

    u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                     un[1:-1, 1:-1] * dt / dx *
                     (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy *
                     (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                     dt / (2 * rho * dx) *
                     (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                     nu * (dt / dx ** 2 *
                           (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                           dt / dy ** 2 *
                           (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])) +
                     F * dt)

    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                     un[1:-1, 1:-1] * dt / dx *
                     (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy *
                     (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                     dt / (2 * rho * dy) *
                     (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                     nu * (dt / dx ** 2 *
                           (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                           dt / dy ** 2 *
                           (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

    # Periodic BC u @ x = 2
    u[1:-1, -1] = (un[1:-1, -1] - un[1:-1, -1] * dt / dx *
                   (un[1:-1, -1] - un[1:-1, -2]) -
                   vn[1:-1, -1] * dt / dy *
                   (un[1:-1, -1] - un[0:-2, -1]) -
                   dt / (2 * rho * dx) *
                   (p[1:-1, 0] - p[1:-1, -2]) +
                   nu * (dt / dx ** 2 *
                         (un[1:-1, 0] - 2 * un[1:-1, -1] + un[1:-1, -2]) +
                         dt / dy ** 2 *
                         (un[2:, -1] - 2 * un[1:-1, -1] + un[0:-2, -1])) + F * dt)

    # Periodic BC u @ x = 0
    u[1:-1, 0] = (un[1:-1, 0] - un[1:-1, 0] * dt / dx *
                  (un[1:-1, 0] - un[1:-1, -1]) -
                  vn[1:-1, 0] * dt / dy *
                  (un[1:-1, 0] - un[0:-2, 0]) -
                  dt / (2 * rho * dx) *
                  (p[1:-1, 1] - p[1:-1, -1]) +
                  nu * (dt / dx ** 2 *
                        (un[1:-1, 1] - 2 * un[1:-1, 0] + un[1:-1, -1]) +
                        dt / dy ** 2 *
                        (un[2:, 0] - 2 * un[1:-1, 0] + un[0:-2, 0])) + F * dt)

    # Periodic BC v @ x = 2
    v[1:-1, -1] = (vn[1:-1, -1] - un[1:-1, -1] * dt / dx *
                   (vn[1:-1, -1] - vn[1:-1, -2]) -
                   vn[1:-1, -1] * dt / dy *
                   (vn[1:-1, -1] - vn[0:-2, -1]) -
                   dt / (2 * rho * dy) *
                   (p[2:, -1] - p[0:-2, -1]) +
                   nu * (dt / dx ** 2 *
                         (vn[1:-1, 0] - 2 * vn[1:-1, -1] + vn[1:-1, -2]) +
                         dt / dy ** 2 *
                         (vn[2:, -1] - 2 * vn[1:-1, -1] + vn[0:-2, -1])))

    # Periodic BC v @ x = 0
    v[1:-1, 0] = (vn[1:-1, 0] - un[1:-1, 0] * dt / dx *
                  (vn[1:-1, 0] - vn[1:-1, -1]) -
                  vn[1:-1, 0] * dt / dy *
                  (vn[1:-1, 0] - vn[0:-2, 0]) -
                  dt / (2 * rho * dy) *
                  (p[2:, 0] - p[0:-2, 0]) +
                  nu * (dt / dx ** 2 *
                        (vn[1:-1, 1] - 2 * vn[1:-1, 0] + vn[1:-1, -1]) +
                        dt / dy ** 2 *
                        (vn[2:, 0] - 2 * vn[1:-1, 0] + vn[0:-2, 0])))

    # Wall BC: u,v = 0 @ y = 0,2
    u[0, :] = 0; v[0, :] = 0
    u[-1, :] = 0; v[-1, :] = 0
    placement = int(nx/2)
    elevation = int(ny/1.4)
    triangleStart = int(nx/3)
    base = int(ny/2)
    triangleHeight = elevation - base
    u[base:triangleHeight, placement] = 0
    u[base, triangleStart:placement] = 0
    tilt = triangleHeight / (placement - triangleStart)

    for it in range(placement - triangleStart + 1):
        u[base:int(base + it * tilt + 1), triangleStart + it] = 0
        u[(int(base - it * tilt - 1)):base, triangleStart + it] = 0
        v[base:int(base + it * tilt + 1), triangleStart + it] = 0
        v[(int(base - it * tilt - 1)):base, triangleStart + it] = 0


    udiff = (np.sum(u) - np.sum(un)) / np.sum(u)
    stepcount += 1
    print(stepcount)


plt.figure(figsize=(24, 4))
X, Y = np.meshgrid(x, y)

plt.plot([triangleStart *dx, placement*dx], [base*dy, (elevation-1)*dy], color='grey', linewidth=2.0)
plt.plot([placement*dx, placement*dx], [base*dy-triangleHeight*dy, (elevation-1)*dy], color='grey', linewidth=2.0)
plt.plot([triangleStart *dx, placement*dx],[base*dy, base*dy-triangleHeight*dy], color='grey', linewidth=2.0)

plt.streamplot(X, Y, u, v, linewidth=1.75, color=np.sqrt(u**2+v**2), cmap=plt.cm.jet, density=2.5)
colorbar = plt.colorbar()
colorbar.set_label('velocity')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

plt.figure(figsize=(24, 4))
X, Y = np.meshgrid(x, y)
plt.contourf(X, Y, np.sqrt(v**2+u**2),150, cmap='jet')
colorbar = plt.colorbar()
plt.plot([triangleStart *dx, placement*dx], [base*dy, (elevation-1)*dy], color='grey', linewidth=2.0)
plt.plot([placement*dx, placement*dx], [base*dy-triangleHeight*dy, (elevation-1)*dy], color='grey', linewidth=2.0)
plt.plot([triangleStart *dx, placement*dx],[base*dy, base*dy-triangleHeight*dy], color='grey', linewidth=2.0)
colorbar.set_label('velocity')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
'''
plt.figure(figsize=(24, 4))
X, Y = np.meshgrid(x, y)
plt.contourf(X, Y, p, cmap='jet')
plt.plot([triangleStart *dx, placement*dx], [base*dy, (elevation-1)*dy], color='grey', linewidth=2.0)
plt.plot([placement*dx, placement*dx], [base*dy-triangleHeight*dy, (elevation-1)*dy], color='grey', linewidth=2.0)
plt.plot([triangleStart *dx, placement*dx],[base*dy, base*dy-triangleHeight*dy], color='grey', linewidth=2.0)
colorbar = plt.colorbar()
colorbar.set_label('p')
plt.streamplot(X, Y, u, v, linewidth=1.75, color='black')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
'''
