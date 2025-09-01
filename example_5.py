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

traj_piecewise = Piecewise((1.54677, x <= 60), (1.63071 - 0.1*x, x > 60))
print("Plotting trajectory. Close plot to continue example...\n")
plot(traj_piecewise, (x, 0, 120))

domain = Interval(48, 69)
xbounds = [45, 72]
ybounds = [-6, 6]
example_name = "Sine-Linear with Square"

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

plot_condition(x, y, cond, xbounds, ybounds, title=example_name)
mathematica_output = print_mathematica(
    x, y, cond, xbounds, ybounds, traj_piecewise, square
)
print("Mathematica command for plotting:\n", mathematica_output)
