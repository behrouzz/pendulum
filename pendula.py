import numpy as np
from scipy.integrate import solve_ivp #odeint


def simple_num(l, th0, vx0=0, vy0=0, dt=0.01, g=9.81, steps=1000):
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



def simple(l, th0, thp0=0, dt=0.01, g=9.81, steps=5000):
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


def elastic(l, k, m, th0=np.pi/2, thp0=0, x0=0, xp0=0, dt=0.01, g=9.81, steps=5000):
    def fun(t, s):
        th = s[0]
        thp = s[1]
        x = s[2]
        xp = s[3]
        thpp = (-g*np.sin(th) - 2*thp*xp) / (l+x)
        xpp = (l+x)*thp**2 - (k/m)*x + g*np.cos(th)
        return thp, thpp, xp, xpp

    t = np.linspace(0, steps*dt, steps)
    t_span = [t[0], t[-1]]
    y0 = [th0, thp0, x0, xp0]
    
    sol = solve_ivp(fun=fun,
                    t_span=t_span,
                    y0=y0,
                    t_eval=t)

    theta = sol.y[0]
    theta_p = sol.y[1]
    x = sol.y[2]
    x_p = sol.y[3]
    
    xx = (l+x) * np.sin(theta)
    yy = (l+x) * np.cos(theta)
    return xx, yy


def elastic_num(l, k, m, x0=1, y0=0, vx0=0, vy0=0, dt=0.01, g=9.81, steps=5000):
    x, y, vx, vy = x0, y0, vx0, vy0

    xx = [x]
    yy = [y]

    for i in range(steps):
        s = np.sqrt(x**2 + y**2) - l
        th = np.arctan(x/y)
        Fx = -k*s*np.sin(th)
        Fy = m*g - k*s*np.cos(th)

        vx = vx + (Fx/m)*dt
        x = x + vx*dt
        vy = vy + (Fy/m)*dt
        y = y + vy*dt

        xx.append(x)
        yy.append(y)

    return np.array(xx), np.array(yy)
