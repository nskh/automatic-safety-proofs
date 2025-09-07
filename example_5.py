#!/usr/bin/env python
# coding: utf-8
from sympy import *
from safe_region_utils import *

init_printing()


x, y = symbols("x y")
## Short Example ##
w = 0.1
square_points: list = [
    geometry.Point(val) for val in [[w, -w], [w, w], [-w, w], [-w, -w]]
]
square: geometry.Polygon = Polygon(*square_points)
print("Plotting object. Close plot to continue example...\n")
plot_polygon(square)

# traj_piecewise = Piecewise((1.54677, x <= 60), (1.63071 - 0.1*x, x > 60))

traj_piecewise = Piecewise((-2.38529, ((x >= 2.0) | (x < 1.0) | (x < 2.0)) & ((x < 1.0) | (x < 2.0) | (x < 3.0)) & ((x >= 1.0) | (x >= 2.0) | (x >= 0)) & ((x >= 1.0) | (x < 3.0) | (x >= 0)) & ((x >= 2.0) | (x < 2.0) | (x >= 0)) & ((x < 2.0) | (x < 3.0) | (x >= 0))), 
          (0.121269*x - 2.7491, (x >= 3.0) & (x < 26.0)),
        (0.0464082*x - 0.802721, (x >= 26.0) & (x <= 27.0)))

        #   (0.121269*x - 2.7491, (x >= 4.0) & (x < 5.0)), 
        #   (0.121269*x - 2.7491, ((x >= 5.0) | (x >= 6.0)) & ((x >= 5.0) | (x < 7.0)) & ((x < 6.0) | (x < 7.0))),
        #   (0.121269*x - 2.7491, (x >= 7.0) & (x < 8.0)),
        # (0.121269*x - 2.7491, ((x >= 8.0) | (x >= 9.0)) & ((x >= 8.0) | (x < 10.0)) & ((x < 9.0) | (x < 10.0))), 
        # (0.121269*x - 2.7491, (x >= 10.0) & (x < 11.0)), (0.121269*x - 2.7491, (x >= 11.0) & (x < 12.0)), 
        # (0.121269*x - 2.7491, (x >= 12.0) & (x < 13.0)), (0.121269*x - 2.7491, (x >= 13.0) & (x <= 14.0)),
        # (0.121269*x - 2.7491, (x >= 14.0) & (x < 15.0)), (0.121269*x - 2.7491, ((x >= 15.0) | (x >= 16.0)) & ((x >= 15.0) | (x < 17.0)) & ((x < 16.0) | (x < 17.0))), 
        # (0.121269*x - 2.7491, (x >= 17.0) & (x < 18.0)), (0.121269*x - 2.7491, (x >= 18.0) & (x < 19.0)), 
        # (0.121269*x - 2.7491, (x >= 19.0) & (x < 20.0)), (0.121269*x - 2.7491, (x >= 20.0) & (x < 21.0)), 
        # (0.121269*x - 2.7491, (x >= 21.0) & (x < 22.0)), (0.121269*x - 2.7491, (x >= 22.0) & (x < 23.0)), 
        # (0.121269*x - 2.7491, (x >= 23.0) & (x < 24.0)), (0.121269*x - 2.7491, ((x >= 24.0) | (x >= 25.0)) & ((x >= 24.0) | (x < 26.0)) & ((x < 25.0) | (x < 26.0))), 
print("Plotting trajectory. Close plot to continue example...\n")
plot(traj_piecewise, (x, 0, 27))

domain = Interval(-8,30)
xbounds = [0, 27]
ybounds = [-3, 2]
example_name = "Relative Trajectory of Two Aircrafts"

cond = compute_unsafe_cond(
    x,
    y,
    square,
    traj_piecewise,
    domain,
)
print("Boolean condition for unsafe region:\n", cond)

print(
    "\nPlotting dot grid visualization of safe and unsafe regions. This may take up to 20 seconds to plot.\nOnce displayed, close plot to continue example...\n"
)

plot_condition(x, y, cond, xbounds, ybounds, title=example_name, resolution=0.1)
mathematica_output = print_mathematica(
    x, y, cond, xbounds, ybounds, traj_piecewise, square
)
print("Mathematica command for plotting:\n", mathematica_output)
