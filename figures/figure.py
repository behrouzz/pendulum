import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

r = 1
theta = np.pi/6


theta1=200
theta2=340

arc = Arc((0, 0), 2*r, 2*r, angle=0,
          theta1=240, theta2=340, fill=False, alpha=0.5)

xp = r * np.sin(theta)
yp = -r * np.cos(theta)


# tangent line at point p
x_tan = np.linspace(0,1)
y_tan = (r**2-xp*x_tan)/yp

# line along spring
m = yp / xp # slope
x_long = np.linspace(xp, xp+0.3)
y_long = m * x_long

fig, ax = plt.subplots()
ax.add_patch(arc)

ax.scatter([xp], [yp], s=100, c='r', marker='s', zorder=10)
ax.plot([-0.1, 0.1], [0, 0], lw=5, c='k')
ax.plot([0, xp], [0, yp], lw=3, c='b')

ax.plot(x_tan, y_tan, lw=1, ls='--', c='r')
ax.plot(x_long, y_long, lw=1, ls='--', c='r')

ax.arrow(xp, yp, dx=0, dy=-0.3, head_width=0.05, length_includes_head=True,
         #head_length=0.1,
         fc='gray',
         ec='gray',
         lw=5, alpha=0.7)

ax.arrow(xp, yp,
         dx=-0.15, dy=-0.15*m, head_width=0.05, length_includes_head=True,
         fc='gray',
         ec='gray',
         lw=5, alpha=0.7)

##ax.annotate('R', (10, 15.6), fontsize=12, weight='bold')
ax.vlines(x=0, ymin=-1, ymax=0, ls='--', lw=1, color='k', zorder=10)
ax.set_aspect('equal')
ax.set_xlim(-0.75, 1.25)
ax.set_ylim(-r-0.5, 0.25)
#plt.axis('off')
#plt.grid()
plt.show()
