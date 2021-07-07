from utils import *
if __name__ == '__main__':
    r = 1
    hexagon = RegularPolygon(Point(0, 0), r, 6)
    # plot_polygon(hexagon)
    diamond = RegularPolygon(Point(0, 0), r, 4)
    # plot_polygon(diamond)
    w = 0.5
    square_points: list = [geometry.Point(val) \
                           for val in [[w, -w], [w, w], [-w, w], [-w, -w]]]
    square: geometry.Polygon = Polygon(*square_points)
    # plot_polygon(square)
    rect_points: list = [geometry.Point(val) \
                         for val in [[2 * w, -w], [2 * w, w], [-2 * w, w], [-2 * w, -w]]]
    rectangle: geometry.Polygon = Polygon(*rect_points)
    # plot_polygon(rectangle)

    x, y = symbols("x y")
    traj_r = 4
    circle_traj = x ** 2 + y ** 2 - traj_r ** 2  # radius 4
    sin_traj = sin(x) - y

    print(safe(square, circle_traj, (4,0), (4.5, 0), x, y), False)
    print(safe(square, circle_traj, (4,0), (4, 0), x, y), False)
    print(safe(square, circle_traj, (4, 0), (3, 0), x, y), True)
    print(safe(square, circle_traj, (4, 0), (5, 0), x, y), True)
    print(safe(diamond, circle_traj, (4, 0), (3, 0), x, y), False)
    print(safe(diamond, circle_traj, (4, 0), (4, 0), x, y), False)
    print(safe(diamond, circle_traj, (4, 0), (2, 0), x, y), True)
    # plot_safe_grid(diamond, circle_traj, (-6, 6), (-6, 6), "Circle traje
    # ctory with diamond", savefig = False)
    # plot_safe_grid(square, circle_traj, (-5, 5), (-5, 5), "Circle trajectory with square")
    # plot_safe_grid(rectangle, sin_traj, (-6, 6), (-3, 3), "Sin trajectory with rectangle")
    # plot_safe_grid(hexagon, sin_traj, (-6, 6), (-3, 3), "Sin trajectory with hexagon")
    # print(outside_active_corners(diamond, circle_traj, Point(0,4), Point(0,6), x, y), True)
    # print(outside_active_corners(diamond, circle_traj, Point(0,4), Point(0,4), x, y), False)
    # print(outside_active_corners(diamond, circle_traj, Point(0,4), Point(0,5), x, y), False)