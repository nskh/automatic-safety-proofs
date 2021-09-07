from utils import *
from time import time
from sympy import solveset, S


def compare_point_in_polygon(poly):
    N = 100
    points = np.random.rand(N, 2)

    start_time = time()
    ray_result = [ray_method(poly, Point(0, 0), Point(p[0], p[1])) for p in points]
    print("Ray: " + str(time() - start_time))

    start_time = time()
    encloses_result = [
        encloses_method(poly, Point(0, 0), Point(p[0], p[1])) for p in points
    ]
    print("Encloses: " + str(time() - start_time))

    assert ray_result == encloses_result


def basic_safe_region_test(square, diamond, circle_traj):
    print(safe(square, circle_traj, Point(4, 0), Point(4.5, 0), x, y), False)
    print(safe(square, circle_traj, Point(4, 0), Point(4, 0), x, y), False)
    print(safe(square, circle_traj, Point(4, 0), Point(3, 0), x, y), True)
    print(safe(square, circle_traj, Point(4, 0), Point(5, 0), x, y), True)
    print(safe(diamond, circle_traj, Point(4, 0), Point(3, 0), x, y), False)
    print(safe(diamond, circle_traj, Point(4, 0), Point(4, 0), x, y), False)
    print(safe(diamond, circle_traj, Point(4, 0), Point(2, 0), x, y), True)


def test_transitions(poly, trajectory, domain, x=None, y=None):
    angles, vertex_pairs = compute_polygon_angles(poly)
    print("Angles:")
    pprint(angles)
    dict_of_transitions, set_of_transitions = find_transitions(
        trajectory, angles, x, y, domain=domain
    )
    print("Set of transitions:")
    pprint(set_of_transitions)
    print("Dict of transitions:")
    pprint(dict_of_transitions)


if __name__ == "__main__":
    r = 1
    hexagon = RegularPolygon(Point(0, 0), r, 6)
    # plot_polygon(hexagon)
    diamond = RegularPolygon(Point(0, 0), r, 4)
    # plot_polygon(diamond)
    w = 0.5
    square_points: list = [
        geometry.Point(val) for val in [[w, -w], [w, w], [-w, w], [-w, -w]]
    ]
    square: geometry.Polygon = Polygon(*square_points)
    # plot_polygon(square)
    rect_points: list = [
        geometry.Point(val)
        for val in [[2 * w, -w], [2 * w, w], [-2 * w, w], [-2 * w, -w]]
    ]
    rectangle: geometry.Polygon = Polygon(*rect_points)
    # plot_polygon(rectangle)

    x, y = symbols("x y")
    traj_r = 4
    circle_traj = x ** 2 + y ** 2 - traj_r ** 2  # radius 4
    sin_traj = sin(x) - y

    # compare_point_in_polygon(diamond)
    # basic_safe_region_test(square, diamond, circle_traj)
    test_transitions(square, circle_traj, domain=Interval(-4, 4))
    test_transitions(rectangle, sin_traj, domain=Interval(-6, 6))

    plot_safe_grid(
        diamond,
        circle_traj,
        (-6, 6),
        (-6, 6),
        "Circle trajectory with diamond",
        savefig=False,
    )
    plot_safe_grid(
        square,
        circle_traj,
        (-5, 5),
        (-5, 5),
        "Circle trajectory with square and notch check",
        Interval(-4, 4),
        savefig=True,
    )
    plot_safe_grid(
        rectangle,
        sin_traj,
        (-6, 6),
        (-3, 3),
        "Sin trajectory with rectangle and notch check",
        Interval(-6, 6),
        savefig=True,
    )
    plot_safe_grid(hexagon, sin_traj, (-6, 6), (-3, 3), "Sin trajectory with hexagon")
