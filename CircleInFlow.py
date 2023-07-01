import numpy as np
import matplotlib.pyplot as plt

nx = 200; ny = 200; nt = 5000; nit = 50; rho = 1; nu = .1; F = 2; dt = .0001; udiff = 1;

L = 2; H = 2

dx = L / (nx - 1); dy = H / (ny - 1); x = np.linspace(0, L, nx); y = np.linspace(0, H, ny)
u = np.zeros((ny, nx)); v = np.zeros((ny, nx)); p = np.ones((ny, nx)); pn = np.ones((ny, nx))
b = np.zeros((ny, nx)); xcc = np.zeros((ny, nx)); ycc = np.zeros((ny, nx))

def circle(radius, centerx, centery, nx):
    theta = np.linspace(0, 2*np.pi, nx)
    x = centerx + radius * np.cos(theta)
    y = centery + radius * np.sin(theta)
    return x, y

radius = 0.21
center = (L/2, H/2)
circ=circle(radius, center[0], center[1], nx)
gridI=np.ones(nx)
gridJ=np.ones(ny)
xc = circ[0]
yc = circ[1]

for i in range(nt):
    un = u.copy(); vn = v.copy(); b = np.zeros_like(u)
    b[1:-1, 1:-1] =(rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) /
    (2 * dy)) - ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx)) ** 2 - 2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) /
    (2 * dy) * (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) - ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) ** 2))

    # Periodic BC Pressure @ x = 2
    b[1:-1, -1] = (rho * (1 / dt * ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx) + (v[2:, -1] - v[0:-2, -1]) / (2 * dy)) -
                  ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx)) ** 2 - 2 * ((u[2:, -1] - u[0:-2, -1]) / (2 * dy) *
                  (v[1:-1, 0] - v[1:-1, -2]) / (2 * dx)) - ((v[2:, -1] - v[0:-2, -1]) / (2 * dy)) ** 2))

    # Periodic BC Pressure @ x = 0
    b[1:-1, 0] = (rho * (1 / dt * ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx) +(v[2:, 0] - v[0:-2, 0]) / (2 * dy))-
                         ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx)) ** 2 -2 * ((u[2:, 0] - u[0:-2, 0]) / (2 * dy) *
                         (v[1:-1, 1] - v[1:-1, -1]) / (2 * dx)) - ((v[2:, 0] - v[0:-2, 0]) / (2 * dy)) ** 2))


    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy ** 2 + (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx ** 2) /
                         (2 * (dx ** 2 + dy ** 2)) - dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) * b[1:-1, 1:-1])

        # Periodic BC Pressure @ x = 2
        p[1:-1, -1] = (((pn[1:-1, 0] + pn[1:-1, -2]) * dy ** 2 +(pn[2:, -1] + pn[0:-2, -1]) * dx ** 2) /
                       (2 * (dx ** 2 + dy ** 2)) - dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) * b[1:-1, -1])

        # Periodic BC Pressure @ x = 0
        p[1:-1, 0] = (((pn[1:-1, 1] + pn[1:-1, -1]) * dy ** 2 + (pn[2:, 0] + pn[0:-2, 0]) * dx ** 2) /
                      (2 * (dx ** 2 + dy ** 2)) -dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) * b[1:-1, 0])

        # Wall boundary conditions, pressure
        p[-1, :] = p[-2, :]; p[0, :] = p[1, :]   # dp/dy = 0 at y = 2 # dp/dy = 0 at y = 0


    u[1:-1, 1:-1] = (un[1:-1, 1:-1] -un[1:-1, 1:-1] * dt / dx *(un[1:-1, 1:-1] - un[1:-1, 0:-2]) -vn[1:-1, 1:-1] * dt / dy *
                    (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -dt / (2 * rho * dx) *(p[1:-1, 2:] - p[1:-1, 0:-2]) +
                    nu * (dt / dx ** 2 *(un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                    dt / dy ** 2 *(un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])) +F * dt)

    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy *(vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) - dt / (2 * rho * dy) *
                     (p[2:, 1:-1]-p[0:-2,1:-1])+nu*(dt/dx**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1, 0:-2]) +
                      dt / dy ** 2 * (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

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
    u[1:-1, 0] = (un[1:-1, 0] - un[1:-1, 0] * dt / dx *(un[1:-1, 0] - un[1:-1, -1]) -
                  vn[1:-1, 0] * dt / dy *(un[1:-1, 0] - un[0:-2, 0]) -
                  dt / (2 * rho * dx) *(p[1:-1, 1] - p[1:-1, -1]) +
                  nu * (dt / dx ** 2 *(un[1:-1, 1] - 2 * un[1:-1, 0] + un[1:-1, -1]) +
                  dt / dy ** 2 *(un[2:, 0] - 2 * un[1:-1, 0] + un[0:-2, 0])) + F * dt)

    # Periodic BC v @ x = 2
    v[1:-1, -1] = (vn[1:-1, -1] - un[1:-1, -1] * dt / dx * (vn[1:-1, -1] - vn[1:-1, -2]) -vn[1:-1, -1] * dt / dy *
                   (vn[1:-1, -1] - vn[0:-2, -1]) -dt / (2 * rho * dy) *(p[2:, -1] - p[0:-2, -1]) +nu *
                   (dt / dx ** 2 *(vn[1:-1, 0] - 2 * vn[1:-1, -1] + vn[1:-1, -2]) +
                    dt / dy ** 2 *(vn[2:, -1] - 2 * vn[1:-1, -1] + vn[0:-2, -1])))

    # Periodic BC v @ x = 0
    v[1:-1, 0] = (vn[1:-1, 0] - un[1:-1, 0] * dt / dx * (vn[1:-1, 0] - vn[1:-1, -1]) -
                  vn[1:-1, 0] * dt / dy *(vn[1:-1, 0] - vn[0:-2, 0]) -dt / (2 * rho * dy) *(p[2:, 0] - p[0:-2, 0]) +
                  nu * (dt / dx ** 2 *(vn[1:-1, 1] - 2 * vn[1:-1, 0] + vn[1:-1, -1]) +
                        dt / dy ** 2 *(vn[2:, 0] - 2 * vn[1:-1, 0] + vn[0:-2, 0])))

    # Wall BC: u,v = 0 @ y = 0,2
    u[0, :] = 0; u[-1, :] = 0
    v[0, :] = 0; v[-1, :] = 0
    for it in range(nx):
        gridI[it] = xc[it] / dx
        gridJ[it] = yc[it] / dy
        u[int(gridJ[it]), int(gridI[it])] = 0
        v[int(gridJ[it]), int(gridI[it])] = 0

    for vert in np.linspace(0,np.pi, nx):
        xglob = center[0] + radius*np.cos(vert)
        yglob = center[1] + radius*np.sin(vert)
        u[int(center[1]/dy):int(yglob/dy),int(xglob/dx)]=0
        v[int(center[1] / dy):int(yglob / dy), int(xglob / dx)] = 0
        xglob = center[0] - radius*np.cos(vert)
        yglob = center[1] - radius*np.sin(vert)
        u[int(yglob/dy):int(center[1]/dy),int(xglob/dx)]=0
        v[int(yglob / dy):int(center[1] / dy), int(xglob / dx)] = 0




    udiff = (np.sum(u) - np.sum(un)) / np.sum(u)
    if udiff <=0.001:
        print("Converged")
        break




plt.figure(figsize=(10, 8))
X, Y = np.meshgrid(x, y)
plt.plot(xc, yc, linewidth=2.5, color='grey')
plt.contourf(X, Y, np.sqrt(v**2+u**2),150, cmap='jet')
colorbar = plt.colorbar()
colorbar.set_label('velocity')
#plt.streamplot(X, Y, u, v, linewidth=1.75, color='black')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

fig = plt.figure(figsize=(10, 8))
ax = plt.gca()
ax.set_facecolor('black')
ax.set_facecolor((0, 0.0, 0.0))
X, Y = np.meshgrid(x, y)
plt.plot(xc, yc, linewidth=2.5, color='grey')
plt.streamplot(X, Y, u, v, linewidth=1.75, color=np.sqrt(u**2+v**2), cmap=plt.cm.jet, density=3.0)
colorbar = plt.colorbar()
colorbar.set_label('velocity')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

plt.figure(figsize=(10, 8))
X, Y = np.meshgrid(x, y)
plt.plot(xc, yc, linewidth=2.5, color='grey')
plt.contourf(X, Y, p,150, cmap='jet')
colorbar = plt.colorbar()
colorbar.set_label('Pressure')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()


