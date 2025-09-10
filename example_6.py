#!/usr/bin/env python
# coding: utf-8
from sympy import *
from safe_region_utils import *

init_printing()


x, y = symbols("x y")
## Short Example ##
w = 1000
square_points: list = [
    geometry.Point(val) for val in [[w, -w], [w, w], [-w, w], [-w, -w]]
]
square: geometry.Polygon = Polygon(*square_points)
print("Plotting object. Close plot to continue example...\n")
plot_polygon(square)

# traj_piecewise = Piecewise((1.54677, x <= 60), (1.63071 - 0.1*x, x > 60))

traj_piecewise = Piecewise((0.0589984962598793*x + 12703.024877962667, x <= -147512.65590951982),
                           (4000.0, (x > -147512.65590951982) & (x <= 73941.3580322652)),
                           (0.022578050905477744*x + 2330.5482543273606, (x > 73941.3580322652) & (x <= 118232.1608206222)),
                           (5000.0, (x > 118232.1608206222) & (x <= 826885.0054343628)),
                           (0.022578050905477737*x + -13669.45174567328, (x > 826885.0054343628) & (x <= 871175.8082227198)), 
                           (6000.0, (x > 871175.8082227198) & (x <= 871175.8082227198)))


     
print("Plotting trajectory. Close plot to continue example...\n")
plot(traj_piecewise, (x, -164462.24044799947, 871175.8082227198))

domain = Interval(-170000, 880000)
xbounds = [-164462.24044799947, 871175.8082227198]
ybounds = [2000, 7000]
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

plot_condition(x, y, cond, xbounds, ybounds, title=example_name, resolution=100)
mathematica_output = print_mathematica(
    x, y, cond, xbounds, ybounds, traj_piecewise, square
)
print("Mathematica command for plotting:\n", mathematica_output)
