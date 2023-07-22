import numpy as np
from animation import animate
from pendula import simple, simple_num, elastic, elastic_num

l = 1
k = 50
m = 1
g = 9.81
th0 = np.pi/4
steps = 1000

xx, yy = elastic(l, k, m, th0=th0, thp0=0, x0=0, xp0=0, steps=steps)
#xx, yy = simple(l=l, th0=th0, steps=steps)
#xx, yy = elastic_num(l, k, m, x0=l, y0=l)

animate(xx, yy)
