import numpy as np
from funcs import animate, el_pend_anal

l = 1
k = 50
m = 1
g = 9.81


xx, yy = el_pend_anal(l, k, m, th0=np.pi/2, thp0=0, x0=0, xp0=0)
animate(xx, yy)
