from utils import *

x, y = symbols("x y")

s = 1  # square side length
# square_points: list = [geometry.Point(val) \
#                        for val in [[0,0], [s, 0], [s, s], [0, s]]]
square_points: list = [
    geometry.Point(val) for val in [[s, 0], [0, s], [-s, 0], [0, -s]]
]
square: geometry.Polygon = Polygon(*square_points)

angles, vertex_pairs = compute_polygon_angles(square)
print(angles)

r = 4
circle_traj = x ** 2 + y ** 2 - r ** 2  # radius 4
plot_implicit(circle_traj)

# TODO(nishant): plot shape on trajectory at given points
print(slope_sym(circle_traj, x, y))

transitions = find_transitions(circle_traj, angles, x, y)
print(transitions)
