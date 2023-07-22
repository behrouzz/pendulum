import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def animate(xx, yy, interval=20):
    fig = plt.figure()
    ax = fig.add_subplot(aspect='equal')

    x0, y0 = xx[0], yy[0]
    line, = ax.plot([0, x0], [0, y0], lw=3, c='k')

    r = np.sqrt(np.abs(xx).max()**2 + np.abs(yy).max()**2)

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
