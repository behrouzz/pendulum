import numpy as np
from astropy import units as u
from funcs import animate, el_pend_num

l = 1 * u.m
k = 50 * u.Unit('N/m')
m = 1 * u.kg
g = 9.81 * u.Unit('m/s2')
dt = 0.01 * u.s

# initial
x0 = l
y0 = 0 * u.m
vx0 = 0 * u.Unit('m/s')
vy0 = 0 * u.Unit('m/s')


xx, yy = el_pend_num(l, k, m, x0=x0, y0=y0, vx0=vx0, vy0=vy0, dt=dt)
animate(xx.value, yy.value)
