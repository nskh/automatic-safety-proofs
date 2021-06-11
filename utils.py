import matplotlib.pyplot as plt
import itertools
import numpy as np
import sympy
from sympy import *

init_printing(use_unicode=True)


def plot_polygon(poly: sympy.Polygon):
    # Draw a polygon by plotting vertices as points and edges as lines
    fig = plt.figure()
    ax = fig.gca()

    verts = poly.vertices

    for p in verts:
        ax.scatter(p.x, p.y, c="r")

    for (p, nextp) in zip(verts, verts[1:] + verts[:1]):
        x = np.linspace(float(p.x), float(nextp.x), 100, dtype=float)
        y = np.linspace(float(p.y), float(nextp.y), 100, dtype=float)
        ax.plot(x, y, c="b")

    ax.axis("equal")

    plt.show()


def compute_polygon_angles(poly: sympy.Polygon) -> list:
    # Go around the polygon and compute the angles (relative to horiz axis)
    # of each of its sides
    verts = poly.vertices
    # Stitch together adjacent vertex pairs (wrap around end of list)
    vertex_pairs = list(zip(verts, verts[1:] + verts[:1]))
    # Angles using atan2: always use atan2 to give correct outputs
    angles = [atan2(nextp.y - p.y, nextp.x - p.x) for (p, nextp) in vertex_pairs]

    # Restrict angles to [0, 2pi]
    positive_angles = [angle if angle > 0 else (angle + 2 * pi) for angle in angles]
    assert max(positive_angles) < 2 * pi
    assert min(positive_angles) >= 0
    return positive_angles, vertex_pairs


def eval_slope(traj, point, x, y):
    # Compute the slope of a trajectory *expression* and plug in an (x,y) point
    return slope_sym(traj, x, y).subs(x, point[0]).subs(y, point[1])


def slope_sym(traj, x, y):
    # TODO: figure out sign issue
    # f(x,y) = 0
    # y = sin(x)
    # dy/dx = cos(x)
    #
    # sin(x) - y = 0
    # df/dx = cos(x)
    # df/dy = -1
    # dy/dx = -cos(x)
    # df/dx / df/dy -> dy/dx
    # TODO: does this hold for things that can't be written y = f(x)?
    # tan(-x) = -tan(x)
    slope = -1 * diff(traj, x) / diff(traj, y)
    return slope


def find_transitions(trajectory, angles, x, y) -> dict:
    # NOTE: trajectory is an *expression*, not equation
    transitions = {}
    for angle in angles:
        # Compute slope symbolically
        slope = slope_sym(trajectory, x, y)
        # TODO(nishant): should this be atan2 here?
        # solve(Eq(atan2(dfdy, dfdx), angle))??
        soln = solve(Eq(slope, tan(angle)))
        for elem in soln:
            # Only add if solution exists (real or dict types)
            if type(elem) == dict or elem.is_real:
                if angle in transitions:
                    transitions[angle].append(elem)
                else:
                    transitions[angle] = [elem]
    print(transitions)

    # soln above may still be symbolic, so we have to evaluate the expression

    transition_points = {}
    traj_eqn = Eq(trajectory, 0)
    for angle, solns in transitions.items():
        # TODO: rename pair variable
        for pair in solns:
            # dict if implicit solution - x as f(y)
            if type(pair) == dict:  # x given as f(y)
                # remove x from equation
                eqn_without_x = traj_eqn.subs(pair)
                y_solns = solve(eqn_without_x)

                for y_soln in y_solns:
                    x_soln_eq = Eq(x, pair[x]).subs(y, y_soln)
                    transition_point = Point(x_soln_eq.rhs, y_soln)
                    if angle in transition_points:
                        transition_points[angle].append(transition_point)
                    else:
                        transition_points[angle] = [transition_point]

            # if not a dictionary, we have an exact solution for x
            else:  # exact solution for x, pair is a single element
                y_solns = solve(traj_eqn.subs(x, pair))
                for y_soln in y_solns:
                    transition_point = Point(pair, y_soln)
                    if angle in transition_points:
                        transition_points[angle].append(transition_point)
                    else:
                        transition_points[angle] = [transition_point]

    return transition_points


# traj -> [x(t); y(t)] -> at some T, what is the angle of the tangent to trajectory
# x,y points, may or may not be on the trajectory,

# at any given (x,y) on traj, there's either a notch OR some pair is active,
# and bounds everything else -> suffices to check for notches AND to check
# that the point is outside *all* active-corner pairs

# TODO(nishant): write better in spec
# given obstacle at (x_O, y_O)
# algortihm is check for polygon inclusion at transition points AND check point
# is outside all active-corner-pairs

# TODO(nishant): create github repo

# TODO(elanor): implement "outside all active-corner pairs" test
# assume symmetric polygons and write something to identify the active corner pairs
# try with a hexagon, square, diamond, rectangle


# Implemented the below two functions when I was confused about how to
# identify active corners in between transition points - but we can just
# test safety using all active-corner pairs so this stuff isn't required.


def find_common_corner(angles_to_vertices: dict, angle, direction):
    # I don't think this is actually useful but keeping it around
    assert abs(direction) == 1

    sorted_angles = sorted(angles_to_vertices.keys())

    angle_index = sorted_angles.index(angle)
    next_greatest_angle_index = (angle_index + direction) % len(
        angles
    )  # wrap around 2pi
    next_greatest_angle = sorted_angles[next_greatest_angle_index]

    common_corner: set = set(angles_to_vertices[angle]).intersection(
        set(angles_to_vertices[next_greatest_angle])
    )
    assert len(common_corner) == 1
    common_corner: Point = common_corner.pop()
    return common_corner


def direction_of_traj_angle(trajectory, transition_points, angle, epsilon):
    # I don't think this is actually useful but keeping it around
    assert epsilon != 0, "Epsilon cannot be 0"

    p = transition_points[angle][0]  # TODO: fix hard coding
    m = tan(angle)
    y_on_tangent = p.y + epsilon * m
    y_solns = solve(Eq(trajectory.subs(x, p.x + epsilon), 0))  # two solutions for y

    # find closest to y_on_tangent
    closest_y_soln = min(y_solns, key=lambda v: abs(v - y_on_tangent))

    if epsilon > 0:
        if closest_y_soln > y_on_tangent:
            # theta increases: pick corners accordingly
            direction = +1
        else:
            # theta decreases: pick corners accordingly
            direction = -1
    elif epsilon < 0:
        if closest_y_soln > y_on_tangent:
            # theta decreases: pick corners accordingly
            direction = -1
        else:
            # theta increases: pick corners accordingly
            direction = +1
    return direction
