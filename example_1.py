#!/usr/bin/env python
# coding: utf-8
from sympy import *
from safe_region_utils import *

init_printing()


x, y = symbols("x y")
## Short Example ##
w = 0.5
square_points: list = [
    geometry.Point(val) for val in [[w, -w], [w, w], [-w, w], [-w, -w]]
]
square: geometry.Polygon = Polygon(*square_points)
plot_polygon(square)

traj_piecewise = Piecewise((sin(x / 2), x < 0), (x / 2, x >= 0))
plot(traj_piecewise)

domain = Interval(-12, 9)
xbounds = [-15, 12]
ybounds = [-3, 9]

cond = compute_unsafe_cond(
    x,
    y,
    square,
    traj_piecewise,
    domain,
)
print("Boolean condition for unsafe region:\n", cond)

print("\nPlotting dot grid visualization of safe and unsafe regions...\n")
plot_condition(x, y, cond, xbounds, ybounds)
mathematica_output = print_mathematica(
    x, y, cond, xbounds, ybounds, traj_piecewise, square, False
)
print("Mathematica command for plotting:\n", mathematica_output)

# TODO(nishant): do we need a f(y) example here? not in the paper
