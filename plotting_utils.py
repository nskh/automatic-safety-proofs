import matplotlib.pyplot as plt
import numpy as np
import sympy
from sympy import *


def plot_condition(
    x,
    y,
    cond,
    xbounds=(-5, 5),
    ybounds=(-5, 5),
    resolution=0.25,
    title="Safe Region Plot",
    alpha=1,
    savefig=False,
):
    fig = plt.figure()
    ax = fig.gca()

    nelem = (xbounds[1] - xbounds[0]) * (ybounds[1] - ybounds[0]) / (resolution ** 2)
    # TODO(nishant): progress bar when plotting somehow?
    count = 0
    for x0 in np.arange(xbounds[0], xbounds[1], resolution):
        for y0 in np.arange(ybounds[0], ybounds[1], resolution):
            count += 1
            # TODO(nishant): progress for plotting dots
            # print(f"{count/nelem*100}% \r")
            # intruder = Point(x0, y0)
            is_safe = (~cond).subs([(x, x0), (y, y0)])
            # is_safe = True
            # for (traj1, traj2) in trajs:
            #     if (traj1.subs(x, intruder[0]).subs(y, intruder[1])) * (
            #         traj2.subs(x, intruder[0]).subs(y, intruder[1])
            #     ) <= 0:
            #         is_safe = False
            #         break
            # if not is_safe:
            #     for transition_point in set_of_transitions:
            #         if is_safe and not encloses_method(
            #             poly, transition_point, intruder
            #         ):
            #             is_safe = False
            #             break
            # TODO(nishant): resolution scaling
            if is_safe:
                ax.plot(x0, y0, "bo", alpha=alpha, markersize=resolution * 4)
            else:
                ax.plot(x0, y0, "ro", alpha=alpha, markersize=resolution * 4)

    ax.axis("equal")

    ax.set_title(title)
    if savefig:
        plt.savefig(title)
    plt.show()


def plot_polygon(poly: sympy.Polygon):
    # Draw a polygon by plotting vertices as points and edges as lines
    fig = plt.figure()
    ax = fig.gca()

    verts = poly.vertices

    for p in verts:
        ax.scatter(p.x, p.y, c="r")
        plt.text(p.x + 0.05, p.y + 0.05, f"({p.x},{p.y})")

    for (p, nextp) in zip(verts, verts[1:] + verts[:1]):
        x = np.linspace(float(p.x), float(nextp.x), 100, dtype=float)
        y = np.linspace(float(p.y), float(nextp.y), 100, dtype=float)
        ax.plot(x, y, c="b")

    ax.axis("equal")

    plt.show()


def move_sympyplot_to_axes(p, ax):
    backend = p.backend(p)
    backend.ax = ax
    # Fix for > sympy v1.5
    backend._process_series(backend.parent._series, ax, backend.parent)
    backend.ax.spines["right"].set_color("none")
    backend.ax.spines["bottom"].set_position("zero")
    backend.ax.spines["top"].set_color("none")
    plt.close(backend.fig)
