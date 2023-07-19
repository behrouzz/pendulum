import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy import units as u
from astropy.units.quantity import Quantity as Q
from scipy.integrate import odeint


def el_pend_num(l, k, m, x0=1, y0=0, vx0=0, vy0=0, dt=0.01, g=9.81, steps=5000):
    if not isinstance(l, Q):
        l = l * u.m
    if not isinstance(k, Q):
        k = k * u.Unit('N/m')
    if not isinstance(m, Q):
        m = m * u.kg
    if not isinstance(x0, Q):
        x0 = x0 * u.m
    if not isinstance(y0, Q):
        y0 = y0 * u.m
    if not isinstance(vx0, Q):
        vx0 = vx0 * u.Unit('m/s')
    if not isinstance(vy0, Q):
        vy0 = vy0 * u.Unit('m/s')
    if not isinstance(dt, Q):
        dt = dt * u.s
    if not isinstance(g, Q):
        g = g * u.Unit('m/s2')

    x, y, vx, vy = x0, y0, vx0, vy0

    xx = [x]
    yy = [y]

    for i in range(steps):
        s = np.sqrt(x**2 + y**2) - l
        th = 2*np.pi*u.rad - np.arctan(x/y)
        Fx = -k*s*np.sin(th)
        Fy = -m*g + k*s*np.cos(th)

        vx = vx + (Fx/m)*dt
        x = x + vx*dt
        vy = vy + (Fy/m)*dt
        y = y + vy*dt

        xx.append(x)
        yy.append(y)

    xx = np.array([i.value for i in xx]) * xx[0].unit
    yy = np.array([i.value for i in yy]) * yy[0].unit
    return xx, yy


def el_pend_anal(l, k, m, th0=np.pi/2, thp0=0, x0=0, xp0=0, dt=0.01, g=9.81, steps=5000):
    def model(y, t, l, k, m):
        th, thp, x, xp = y
        thpp = (-g*np.sin(th) - 2*thp*xp) / (l+x)
        xpp = (l+x)*thp**2 - (k/m)*x + g*np.cos(th)
        return thp, thpp, xp, xpp

    t = np.linspace(0, steps*dt, steps)

    y0 = [th0, thp0, x0, xp0]
    y = odeint(model, y0, t, args=(l, k, m))

    theta, x = y[:,0], y[:,2]
    
    xx = (l+x) * np.sin(theta)
    yy = -(l+x) * np.cos(theta)
    return xx, yy


def set_lim(ax, xx, yy):
    minxs = min([xx.min(), 0])
    maxxs = max([xx.max(), 0])
    minys = min([yy.min(), 0])
    maxys = max([yy.max(), 0])

    minxs = minxs*1.2 if minxs<0 else minxs*0.8
    minys = minys*1.2 if minys<0 else minys*0.8
    maxxs = maxxs*0.8 if maxxs<0 else maxxs*1.2
    maxys = maxys*0.8 if maxys<0 else maxys*1.2

    if maxys==0:
        maxys=0.5
    
    ax.set_xlim([minxs, maxxs])
    ax.set_ylim([minys, maxys])
    
    return ax


def animate(xx, yy, interval=20):
    fig = plt.figure()
    ax = fig.add_subplot(aspect='equal')

    x0, y0 = xx[0], yy[0]
    line, = ax.plot([0, x0], [0, y0], lw=3, c='k')

    bob_radius = 0.08
    cir = ax.add_patch(plt.Circle((x0, y0), bob_radius, fc='r', zorder=3))

    ax = set_lim(ax, xx, yy)

    def animate(i):
        x, y = xx[i], yy[i]
        line.set_data([0, x], [0, y])
        cir.set_center((x, y))

    #ax.invert_yaxis()
    ani = animation.FuncAnimation(fig, animate, frames=len(xx), repeat=True, interval=interval)
    plt.show()
