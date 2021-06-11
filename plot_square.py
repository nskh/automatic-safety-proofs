from utils import *

x, y = symbols("x y")

s = 1  # square side length
# square_points: list = [geometry.Point(val) \
#                        for val in [[0,0], [s, 0], [s, s], [0, s]]]
square_points: list = [
    geometry.Point(val) for val in [[s, 0], [0, s], [-s, 0], [0, -s]]
]
square: geometry.Polygon = Polygon(*square_points)
plot_polygon(square)

angles, vertex_pairs = compute_polygon_angles(square)
print(angles)
