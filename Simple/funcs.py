import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy import units as u
from astropy.units.quantity import Quantity as Q
from scipy.integrate import odeint, solve_ivp


def simple_pend_num(l, th0, vx0=0, vy0=0, dt=0.01, g=9.81, steps=1000):
    th, vx, vy = th0*(np.pi/180), vx0, vy0

    x = l * np.sin(th)
    y = l * np.cos(th)
    
    m = 1
    xx = [x]
    yy = [y]

    for i in range(steps):
        th = np.arcsin(x/l)
        Fx = -m*g*np.sin(th)*np.cos(th)
        vx = vx + (Fx/m)*dt
        x = x + vx*dt
        y = np.sqrt(l**2-x**2)
        xx.append(x)
        yy.append(y)

    return np.array(xx), np.array(yy)


##def simple_pend_anal(l, th0, thp0=0, dt=0.01, g=9.81, steps=5000):
##    def model(y, t, l):
##        th, thp = y
##        thpp = -(g/l)*np.sin(th)
##        return thp, thpp
##
##    t = np.linspace(0, steps*dt, steps)
##
##    y0 = [np.radians(th0), thp0]
##    y = odeint(model, y0, t, args=(l,))
##
##    theta = y[:,0]
##    xx = l * np.sin(theta)
##    yy = l * np.cos(theta)
##    return xx, yy



def simple_pend_anal(l, th0, thp0=0, dt=0.01, g=9.81, steps=5000):
    def fun(t, s):
        th = s[0]
        thp = s[1]
        thpp = -(g/l)*np.sin(s[0])
        return thp, thpp

    t = np.arange(0, 10.01, 0.01)
    t_span = [t[0], t[-1]]
    y0 = [th0, thp0] 

    sol = solve_ivp(fun=fun,
                    t_span=t_span,
                    y0=y0,
                    t_eval=t)

    theta = sol.y[0]
    #theta_p = sol.y[1]
    xx = l * np.sin(theta)
    yy = l * np.cos(theta)
    return xx, yy


def animate(xx, yy, interval=20):
    fig = plt.figure()
    ax = fig.add_subplot(aspect='equal')

    x0, y0 = xx[0], yy[0]
    line, = ax.plot([0, x0], [0, y0], lw=3, c='k')

    r = np.sqrt(xx[0]**2 + yy[0]**2)

    radius = 0.05 * r
    cir = ax.add_patch(plt.Circle((x0, y0), radius, fc='r', zorder=3))
    
    ax.set_xlim([-1.2*r, 1.2*r])
    ax.set_ylim([-1.2*r, 1.2*r])

    def anim(i):
        x, y = xx[i], yy[i]
        line.set_data([0, x], [0, y])
        cir.set_center((x, y))

    ax.invert_yaxis()
    ani = animation.FuncAnimation(fig, anim, frames=len(xx), repeat=True, interval=interval)
    plt.show()
