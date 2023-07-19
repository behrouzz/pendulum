import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

r = 1
theta = np.pi/6


theta1=200
theta2=340

arc = Arc((0, 0), 2*r, 2*r, angle=0, theta1=200, theta2=340, fill=False, alpha=0.3)

xp = r * np.sin(theta)
yp = -r * np.cos(theta)

fig, ax = plt.subplots()
ax.add_patch(arc)

ax.scatter([0], [0], s=25, c='k', marker='P')
ax.scatter([xp], [yp], s=50, c='r', marker='o')
ax.plot([-0.1, 0.1], [0, 0], lw=5, c='k')
ax.plot([0, xp], [0, yp])

##ax.annotate('R', (10, 15.6), fontsize=12, weight='bold')
#ax.axvline(x=0, ymin=1.2, ymax=0, ls='--', lw=1, c='k', zorder=10)
ax.vlines(x=0, ymin=-1, ymax=0, ls='--', lw=1, color='k', zorder=10)
ax.set_aspect('equal')
ax.set_xlim(-1, 1)
ax.set_ylim(-r-0.25, 0.25)
plt.axis('off')
#plt.grid()
plt.show()
