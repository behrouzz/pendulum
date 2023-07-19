import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy import units as u
from astropy.units.quantity import Quantity as Q
from scipy.integrate import odeint

def where(th):
    q1 = 0*u.rad           <= th < (np.pi/2)*u.rad
    q2 = (np.pi/2)*u.rad   <= th < np.pi*u.rad
    q3 = np.pi*u.rad       <= th < (3*np.pi/2)*u.rad
    q4 = (3*np.pi/2)*u.rad <= th < 2*np.pi2*u.rad
    return q1, q2, q3, q4


def simple_pend_num(l, th0=45, vx0=0, vy0=0, dt=0.01, g=9.81, steps=1000):
    if not isinstance(l, Q):
        l = l * u.m
    if not isinstance(th0, Q):
        th0 = (th0 * u.deg).to('rad')
    if not isinstance(vx0, Q):
        vx0 = vx0 * u.Unit('m/s')
    if not isinstance(vy0, Q):
        vy0 = vy0 * u.Unit('m/s')
    if not isinstance(dt, Q):
        dt = dt * u.s
    if not isinstance(g, Q):
        g = g * u.Unit('m/s2')

    th, vx, vy = th0, vx0, vy0

    q1, q2, q3, q4 = where(th)

    if q1:
        x = l * np.sin(th)
        y = l * np.cos(th)
    if q2:
        x = l * np.cos(th-(np.pi/2)*u.rad)
        y = -l * np.sin(th-(np.pi/2)*u.rad)
    else:
        print('*'*70)


    m = 1

    xx = [x]
    yy = [y]

    for i in range(steps):    
        print(th.to('deg'), '=>', x)
        if q1:
            th = np.arcsin(x/l)
        elif q2:
            th = np.arccos(x/l)+(np.pi/2)*u.rad
        Fx = -m*g*np.sin(th)*np.cos(th)
        vx = vx + (Fx/m)*dt
        x = x + vx*dt
        y = np.sqrt(l**2-x**2) if x<l else 0
        if q2:
            y = -y
        
        xx.append(x)
        yy.append(y)
        q1, q2, q3, q4 = where(th)

    xx = np.array([i.value for i in xx]) * xx[0].unit
    yy = np.array([i.value for i in yy]) * yy[0].unit
    return xx, yy


def simple_pend_anal(l, th0, thp0=0, dt=0.01, g=9.81, steps=5000):
    def model(y, t, l):
        th, thp = y
        thpp = -(g/l)*np.sin(th)
        return thp, thpp

    t = np.linspace(0, steps*dt, steps)

    y0 = [np.radians(th0), thp0]
    y = odeint(model, y0, t, args=(l,))

    theta = y[:,0]
    xx = l * np.sin(theta)
    yy = l * np.cos(theta)
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

    r = np.sqrt(xx[0]**2 + yy[0]**2)

    radius = 0.05 * r
    cir = ax.add_patch(plt.Circle((x0, y0), radius, fc='r', zorder=3))

    #ax = set_lim(ax, xx, yy)
    
    ax.set_xlim([-1.2*r, 1.2*r])
    ax.set_ylim([-1.2*r, 1.2*r])

    def anim(i):
        x, y = xx[i], yy[i]
        line.set_data([0, x], [0, y])
        cir.set_center((x, y))

    ax.invert_yaxis()
    ani = animation.FuncAnimation(fig, anim, frames=len(xx), repeat=True, interval=interval)
    plt.show()
