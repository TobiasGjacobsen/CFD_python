import numpy as np
import matplotlib.pyplot as plt


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


def pressure_poisson_periodic(p, dx, dy, cornerx, cornery):

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

        p[41:cornery,cornerx]=p[41:cornery, cornerx-1]
        p[cornery, 41:cornerx]=p[cornery-1, 41:cornerx]
        p[41:cornery,39 ] = p[41:cornery, 40]
        p[39, 41:cornerx]=p[40, 41:cornerx]

        '''
        u[41:60,60] = 0; v[41:60,60] = 0

        u[60, 41:60]= 0; v[60, 41:60]= 0;

        u[41:60, 41] = 0; v[41:60, 41] = 0;

        u[41, 41:60] = 0; v[41, 41:60] = 0;
        '''



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
F = 14
dt = .001
cornerx = 80
cornery = 80

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
    p = pressure_poisson_periodic(p, dx, dy, cornerx, cornery)

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
    u[-1, :] = 0;  v[-1, :] = 0


    u[41:cornery,60] = 0; v[41:cornery,60] = 0

    u[60, 41:cornerx]= 0; v[60, 41:cornerx]= 0;

    u[41:cornery, 41] = 0; v[41:cornery, 41] = 0;

    u[41, 41:cornery] = 0; v[41, 41:cornerx] = 0;
    u[41:cornery, 41:cornerx]=0; v[41:cornery, 41:cornerx]=0





    udiff = (np.sum(u) - np.sum(un)) / np.sum(u)
    stepcount += 1
    print(stepcount)



plt.figure(figsize=(24, 4))
X, Y = np.meshgrid(x, y)
plt.contourf(X, Y, np.sqrt(v**2+u**2),150, cmap='jet')
colorbar = plt.colorbar()
colorbar.set_label('velocity')
plt.plot([41*dx,41*dx], [41*dy, cornery*dy], linewidth=3.0, color='grey')
plt.plot([41*dx,cornerx*dx], [cornery*dy, cornery*dy], linewidth=3.0, color='grey')
plt.plot([cornerx*dx,cornerx*dx], [41*dy, cornery*dy], linewidth=3.0, color='grey')
plt.plot([41*dx,cornerx*dx], [41*dy, 41*dy], linewidth=3.0, color='grey')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

plt.figure(figsize=(24, 4))
X, Y = np.meshgrid(x, y)
plt.contourf(X, Y, np.sqrt(v**2+u**2), 150, cmap='jet')
colorbar = plt.colorbar()
colorbar.set_label('velocity')
plt.streamplot(X, Y, u, v, linewidth=1.2, color='black', density=1.75)
plt.plot([41*dx,41*dx], [41*dy, (cornery)*dy], linewidth=3.0, color='grey')
plt.plot([41*dx,(cornerx-1)*dx], [(cornery)*dy, (cornery)*dy], linewidth=2.0, color='grey')
plt.plot([(cornerx-1)*dx,(cornerx-1)*dx], [41*dy, (cornery)*dy], linewidth=2.0, color='grey')
plt.plot([41*dx,(cornerx-1)*dx], [41*dy, 41*dy], linewidth=3.0, color='grey')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

plt.figure(figsize=(24, 4))
X, Y = np.meshgrid(x, y)
plt.streamplot(X, Y, u, v, linewidth=1.75, color=np.sqrt(u**2+v**2), cmap=plt.cm.jet)
colorbar = plt.colorbar()
colorbar.set_label('velocity')
plt.plot([41*dx,41*dx], [41*dy, (cornery)*dy], linewidth=3.0, color='grey')
plt.plot([41*dx,(cornerx-1)*dx], [(cornery)*dy, (cornery)*dy], linewidth=2.0, color='grey')
plt.plot([(cornerx-1)*dx,(cornerx-1)*dx], [41*dy, (cornery)*dy], linewidth=2.0, color='grey')
plt.plot([41*dx,(cornerx-1)*dx], [41*dy, 41*dy], linewidth=3.0, color='grey')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

plt.figure(figsize=(24, 4))
X, Y = np.meshgrid(x, y)
plt.contourf(X, Y, p,150, cmap='jet')
colorbar = plt.colorbar()
colorbar.set_label('p')
plt.plot([41*dx,41*dx], [41*dy, cornery*dy], linewidth=3.0, color='grey')
plt.plot([41*dx,cornerx*dx], [cornery*dy, cornery*dy], linewidth=3.0, color='grey')
plt.plot([cornerx*dx,cornerx*dx], [41*dy, cornery*dy], linewidth=3.0, color='grey')
plt.plot([41*dx,cornerx*dx], [41*dy, 41*dy], linewidth=3.0, color='grey')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

