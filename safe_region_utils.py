import matplotlib.pyplot as plt
import numpy as np
import sympy
from sympy import *
from typing import Tuple, List, Dict, Set

from mathematica_utils import *
from plotting_utils import *

init_printing(use_unicode=True)


def compute_polygon_angles(poly: sympy.Polygon) -> list:
    # Go around the polygon and compute the angles (relative to horiz axis)
    # of each of its sides
    verts = poly.vertices
    # Stitch together adjacent vertex pairs (wrap around end of list)
    vertex_pairs = list(zip(verts, verts[1:] + verts[:1]))
    # Angles using atan2: always use atan2 to give correct outputs
    vertex_offsets: List[Tuple] = [
        (nextp.x - p.x, nextp.y - p.y) for (p, nextp) in vertex_pairs
    ]
    angles: List[float] = [
        atan2(offset_y, offset_x) for (offset_x, offset_y) in vertex_offsets
    ]

    # Restrict angles to [0, 2pi]
    positive_angles = [angle % (2 * pi) for angle in angles]
    # assert max(positive_angles) < 2 * pi
    # assert min(positive_angles) >= 0
    return positive_angles, vertex_pairs


def eval_slope(traj, point, x, y):
    # Compute the slope of a trajectory *expression* and plug in an (x,y) point
    df_dy, df_dx = slope_sym(traj, x, y)
    return df_dy.subs(x, point[0]).subs(y, point[1]), df_dx.subs(x, point[0]).subs(
        y, point[1]
    )


# TODO(nishant): fix, dear god
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
    # slope = -1 * diff(traj, x) / diff(traj, y)
    return -1 * diff(traj, x), diff(traj, y)


def find_transitions(trajectory, angles, x, y, domain=S.Complexes):
    """
    TODO(nishant): add docstrings

    Returns
    =======
    transition_points: dict
    set_of_transitions: set
    """
    # NOTE: trajectory is an *expression*, not equation
    transitions = {}
    # TODO(nishant): explain slope and variables for non f(x,y) functions
    # TODO(nishant): see if we actually need negative df/dy?
    df_dy, df_dx = slope_sym(trajectory, x, y)
    slope = df_dy / df_dx
    for angle in angles:
        # Compute slope symbolically
        # soln = solve(Eq(df_dx * sin(angle), df_dy * cos(angle)), dict=True)
        # 2 vectors <x1, y1> <x2, y2>
        # parallel iff y1*x2 = x1*y2
        # <dfdx, dfdy> <cos(theta), sin(theta)>
        # TODO(nishant): label data types better
        # TODO(nishant): rename variables
        soln = solveset(Eq(df_dx * sin(angle), df_dy * cos(angle)), x, domain=domain)

        if soln == S.EmptySet:
            soln = solveset(
                Eq(df_dx * sin(angle), df_dy * cos(angle)), y, domain=domain
            )
            if soln != S.EmptySet and type(soln) != list:
                # In this case, type(soln): S.FiniteSet
                soln = [{y: soln_elem} for soln_elem in list(soln)]
        else:
            # Pack into list of dict so it's clear which variable has been solved for
            soln = [{x: soln_elem} for soln_elem in list(soln)]

        for elem in soln:
            if angle in transitions:
                transitions[angle].append(elem)
            else:
                transitions[angle] = [elem]
    # soln above may still be symbolic, so we have to evaluate the expression
    # that's what happens below

    transition_points = {}
    set_of_transitions = set()
    traj_eqn = Eq(trajectory, 0)
    for angle, solns in transitions.items():
        for pair in solns:
            # pair should always be a dictionary
            assert type(pair) == dict, "Solution element was not a dictionary!"
            # pair looks like {x: f(y)} or {y: f(x)}
            # remove one variable from equation by substituting pair into traj_eqn
            traj_eqn_single_var = traj_eqn.subs(pair)
            # traj_eqn used to have two variables but now has only one
            single_var_solns = solve(traj_eqn_single_var)

            # before going further, figure out the variable for
            # which pair contains a solution
            soln_var = [k for k in pair][0]  # variable is the dict key
            other_var = y if soln_var == x else x

            for single_var_soln in single_var_solns:
                # substitute in single_var_soln to solve for soln_var
                solved_eqn = Eq(soln_var, pair[soln_var]).subs(
                    other_var, single_var_soln
                )
                # with this, we have a solution for the transition point
                if soln_var == x:
                    transition_point = Point(solved_eqn.rhs, single_var_soln)
                elif soln_var == y:
                    transition_point = Point(single_var_soln, solved_eqn.rhs)
                set_of_transitions.add(transition_point)
                if angle in transition_points:
                    transition_points[angle].append(transition_point)
                else:
                    transition_points[angle] = [transition_point]

    return transition_points, set_of_transitions


def ray_method(poly: sympy.Polygon, location: sympy.Point, intruder: sympy.Point):
    shifted_intruder = Point(intruder) - Point(location)
    intersections = poly.intersection(
        Segment2D(
            shifted_intruder,
            Point(shifted_intruder[0] + poly.perimeter, shifted_intruder[1]),
        )
    )
    return len(intersections) % 2 == 0


def encloses_method(poly: sympy.Polygon, location: sympy.Point, intruder: sympy.Point):
    shifted_intruder: Point = Point(intruder) - Point(location)
    if shifted_intruder in poly.vertices or any(
        shifted_intruder in s for s in poly.sides
    ):
        return False
    return not poly.encloses_point(shifted_intruder)


# TODO(nishant): rename this
def compute_unsafe_cond(
    x,
    y,
    poly: sympy.Polygon,
    trajectory,  # piecewise
    domain,
    add_notches=True,
    print_latex=False,
):
    """Generates clauses
    TODO(nishant): docstring
    """

    # TODO(nishant): figure out input spec
    # TODO(nishant): test piecewise with functions of y

    angles, vertex_pairs = compute_polygon_angles(poly)
    verts = poly.vertices
    # compute width to use for g()
    # TODO(nishant): make this work for Y also
    if y not in trajectory.free_symbols:
        func_var = x
        keyfunc = lambda p: p.x
    elif x not in trajectory.free_symbols:
        func_var = y
        keyfunc = lambda p: p.y
    else:
        raise Exception("Trajectory had two variables!")

    w_point = max(
        [
            -1 * min([v - poly.centroid for v in verts], key=keyfunc),
            max([v - poly.centroid for v in verts], key=keyfunc),
        ],
        key=keyfunc,
    )
    w = getattr(w_point, str(func_var))

    set_of_transitions = set()
    if type(trajectory) == Piecewise:
        for (subtraj, subcond) in trajectory.as_expr_set_pairs():
            # trim domain by computing intersection
            subdomain = subcond.intersect(domain)

            # TODO(nishant) handle f(x) and f(y) trajs

            # TODO(nishant): assert only 1 free symbol
            # find transitions for piecewise domain over this interval
            func_var = subtraj.free_symbols.pop()

            if y not in subtraj.free_symbols:  # form y=f(x)
                _, subset_of_transitions = find_transitions(
                    -y + subtraj, angles, x, y, domain=subdomain
                )
                # add piecewise boundary
                left_bound = Point(subdomain.inf, subtraj.subs(func_var, subdomain.inf))
                right_bound = Point(
                    subdomain.sup, subtraj.subs(func_var, subdomain.sup)
                )
            elif x not in subtraj.free_symbols:  # form x=f(y)
                _, subset_of_transitions = find_transitions(
                    -x + subtraj, angles, x, y, domain=subdomain
                )
                # add piecewise boundary
                left_bound = Point(subtraj.subs(func_var, subdomain.inf), subdomain.inf)
                right_bound = Point(
                    subtraj.subs(func_var, subdomain.sup), subdomain.sup
                )
            set_of_transitions.update(subset_of_transitions)

            set_of_transitions.add(left_bound)
            set_of_transitions.add(right_bound)
    else:
        if y not in trajectory.free_symbols:
            _, subset_of_transitions = find_transitions(
                -y + trajectory, angles, x, y, domain=domain
            )
            left_bound = Point(domain.inf, trajectory.subs(func_var, domain.inf))
            right_bound = Point(domain.sup, trajectory.subs(func_var, domain.sup))
        elif x not in trajectory.free_symbols:
            _, subset_of_transitions = find_transitions(
                -x + trajectory, angles, x, y, domain=domain
            )
            left_bound = Point(trajectory.subs(func_var, domain.inf), domain.inf)
            right_bound = Point(trajectory.subs(func_var, domain.sup), domain.sup)

        set_of_transitions.update(subset_of_transitions)
        set_of_transitions.add(left_bound)
        set_of_transitions.add(right_bound)

    # TODO(nishant): test with f(y)
    sorted_transitions: list = sorted(
        set_of_transitions, key=lambda point: getattr(point, str(func_var))
    )
    func_var_transitions = [getattr(p, str(func_var)) for p in sorted_transitions]
    midpoints = np.convolve(func_var_transitions, [1, 1], "valid") / 2
    if func_var == x:
        # sorted_transitions: list = sorted(set_of_transitions, key=lambda point: point.x)
        # x_transitions = [p.x for p in sorted_transitions]
        # midpoints = np.convolve(x_transitions, [1, 1], "valid") / 2
        deriv_traj = diff(trajectory, x)
        # dydx_midpoints = [deriv_traj.subs(x, val) for val in midpoints]
        # midpoint_angles = [atan2(d, 1) for d in dydx_midpoints]
    elif func_var == y:
        # sorted_transitions: list = sorted(set_of_transitions, key=lambda point: point.y)
        # y_transitions = [p.y for p in sorted_transitions]
        # midpoints = np.convolve(y_transitions, [1, 1], "valid") / 2
        # TODO(nishant): test this for f(y)!
        deriv_traj = 1 / diff(trajectory, y)  # invert dx/dy for slope
        # dydx_midpoints = [deriv_traj.subs(y, val) for val in midpoints]
    dydx_midpoints = [deriv_traj.subs(func_var, val) for val in midpoints]
    midpoint_angles = [atan2(d, 1) for d in dydx_midpoints]

    # TODO(nishant): more comments in this function

    vertex_pairs_right = list(zip(verts, verts[1:] + verts[:1]))
    vertex_offsets_right: List[Tuple] = [
        (nextp.x - p.x, nextp.y - p.y) for (p, nextp) in vertex_pairs_right
    ]
    angles_right: List[float] = [
        atan2(offset_y, offset_x) for (offset_x, offset_y) in vertex_offsets_right
    ]
    positive_angles_right = [angle % (2 * pi) for angle in angles_right]

    vertex_pairs_left = list(zip(verts[-1:] + verts[:-1], verts))
    vertex_offsets_left: List[Tuple] = [
        (nextp.x - p.x, nextp.y - p.y) for (p, nextp) in vertex_pairs_left
    ]
    angles_left: List[float] = [
        atan2(offset_y, offset_x) for (offset_x, offset_y) in vertex_offsets_left
    ]
    positive_angles_left = [angle % (2 * pi) for angle in angles_left]

    corners_to_angles = dict()
    for i in range(len(verts)):
        if positive_angles_right[i] < positive_angles_left[i]:
            corners_to_angles[verts[i]] = Interval(
                positive_angles_left[i], positive_angles_right[i] + 2 * pi
            )
        else:
            corners_to_angles[verts[i]] = Interval(
                positive_angles_left[i], positive_angles_right[i]
            )

    active_corners = {}
    for midpoint_angle in midpoint_angles:
        for k, v in corners_to_angles.items():
            # TODO(nishant): if not symmetric, add 180 deg
            if midpoint_angle % (2 * pi) in v:
                active_corners[midpoint_angle] = k
    var_intervals = list(
        zip(func_var_transitions[:-1], func_var_transitions[1:])
    )  # same len as midpoints

    total_cond = None
    for i, var_interval in enumerate(var_intervals):
        if func_var == x:
            x_left = var_interval[0]
            y_left = trajectory.subs(x, x_left)
            x_right = var_interval[1]
            y_right = trajectory.subs(x, x_right)
        elif func_var == y:
            y_left = var_interval[0]
            x_left = trajectory.subs(y, y_left)
            y_right = var_interval[1]
            x_right = trajectory.subs(y, y_right)

        poly_center_left = Point(x_left, y_left)
        poly_center_right = Point(x_right, y_right)

        # TODO(nishant): if asymmetric use two corners, not center for constructing line
        active_corner_offset = active_corners[midpoint_angles[i]]
        line_left = Line(poly_center_left, poly_center_left + active_corner_offset)
        line_right = Line(poly_center_right, poly_center_right + active_corner_offset)
        # using Line.equation() creates duplicate x,y variables and ruins .subs() call later
        left_a, left_b, left_c = line_left.coefficients
        right_a, right_b, right_c = line_right.coefficients
        left_eq = left_a * x + left_b * y + left_c
        right_eq = right_a * x + right_b * y + right_c
        corner_cond = left_eq * right_eq <= 0

        # TODO(nishant): change this to handle y functions too
        if func_var == x:
            func_var_cond = (x >= (x_left - w)) & (x <= (x_right + w))
        elif func_var == y:
            func_var_cond = (y >= (y_left - w)) & (y <= (y_right + w))

        # construct g functions
        # trajectory over this interval, held constant outside of it
        if type(trajectory) == Piecewise:
            if func_var == x:
                current_piece = trajectory.as_expr_set_pairs(Interval(x_left, x_right))
                # TODO(nishant): may yield two solutions in some situations - left inclusive?
                # need to filter better than just taking 0th element
                current_fn = current_piece[0][0]
                g = y - Piecewise(
                    (y_left, x < x_left),
                    (current_fn, x <= x_right),
                    (y_right, x > x_right),
                )
            elif func_var == y:
                current_piece = trajectory.as_expr_set_pairs(Interval(y_left, y_right))
                # TODO(nishant): may yield two solutions in some situations - left inclusive?
                # need to filter better than just taking 0th element
                current_fn = current_piece[0][0]
                g = x - Piecewise(
                    (x_left, y < y_left),
                    (current_fn, y <= y_right),
                    (x_right, y > y_right),
                )
        else:
            if func_var == x:
                g = y - Piecewise(
                    (y_left, x < x_left),
                    (trajectory, x <= x_right),
                    (y_right, x > x_right),
                )
            elif func_var == y:
                g = x - Piecewise(
                    (x_left, y < y_left),
                    (trajectory, y <= y_right),
                    (x_right, y > y_right),
                )
        # TODO(nishant): rewrite for asymmetric
        g1 = g.subs(x, x - active_corner_offset.x).subs(y, y - active_corner_offset.y)
        g2 = g.subs(x, x + active_corner_offset.x).subs(y, y + active_corner_offset.y)

        full_cond = corner_cond & func_var_cond & (g1 * g2 <= 0)
        if total_cond is None:
            total_cond = full_cond
        else:
            total_cond = total_cond | full_cond

    if add_notches:
        # adding notches
        for transition_point in set_of_transitions:
            # neg for left side, 0 for on edge, pos for on right side
            # inside polygon IF all neg or IF all pos
            shifted_vertex_pairs = [
                (v + transition_point, nextv + transition_point)
                for (v, nextv) in vertex_pairs
            ]  # [(v1, v2), (v2, v3), ..., (vn, v1)]
            # source: https://inginious.org/course/competitive-programming/geometry-pointinconvex
            # source: https://www.eecs.umich.edu/courses/eecs380/HANDOUTS/PROJ2/InsidePoly.html
            side_conds = [
                (y - v.y) * (nextv.x - v.x) - (x - v.x) * (nextv.y - v.y)
                for (v, nextv) in shifted_vertex_pairs
            ]
            # construct both cases inside
            cond1 = true  # init to sympy.true since we're cascading Ands
            cond2 = true  # init to sympy.true since we're cascading Ands
            for side_cond in side_conds:
                cond1 = cond1 & (side_cond <= 0)
                cond2 = cond2 & (side_cond >= 0)
            total_cond = total_cond | cond1 | cond2

    if print_latex:
        print(latex(total_cond))

    return total_cond


# TODO(nishant): delete this function after checking its functionality is duplicated above
def plot_safe_grid(
    poly: sympy.Polygon,
    trajectory,
    xbounds,
    ybounds,
    title,
    domain,
    resolution=0.25,
    alpha=1,
    savefig=False,
    add_notches=True,
    piecewise_bounds=None,
):
    fig = plt.figure()
    ax = fig.gca()
    # p1 = plot_implicit(trajectory, line_color='k')
    # backend = p1.backend(p1)
    # ax = backend.ax

    x, y = symbols("x y")
    verts = poly.vertices
    trajs = []
    # assume even number of vertices and symmetric polygon
    assert len(verts) % 2 == 0, "Polygon given did not have an even number of sides!"
    for i in range(int(len(verts) / 2)):
        # assume center of polygon is (0,0)
        offset = verts[i]  # - poly.center
        trajs.append(
            (
                trajectory.subs(x, x - offset[0]).subs(y, y - offset[1]),
                trajectory.subs(x, x + offset[0]).subs(y, y + offset[1]),
            )
        )
    angles, vertex_pairs = compute_polygon_angles(poly)

    if type(trajectory) is list:
        # Then we're operating with a piecewise function
        assert type(domain) is list, "type(domain) was not a list"
        assert len(trajectory) == len(
            domain
        ), "trajectories and domains were not of the same length"
        set_of_transitions = set()
        for subtraj, subdomain in zip(trajectory, domain):
            _, subset_of_transitions = find_transitions(
                subtraj, angles, x, y, domain=subdomain
            )
            set_of_transitions.update(subset_of_transitions)

        # add each piecewise boundary also
        # TODO(nishant): add better checking this coheres with domains, or automate
        assert piecewise_bounds is not None, "must specify piecewise_bounds!"
        set_of_transitions.update(piecewise_bounds)

        print(set_of_transitions)
    else:
        dict_of_transitions, set_of_transitions = find_transitions(
            trajectory, angles, x, y, domain=domain
        )

    # order each transition point [begin with left bound, end w right bound]
    # in between each boundary point find the midpoint and identify active corners
    # for each boundary:
    #   compute f(x) and store in dict mapping x_bound -> f(x_bound)
    #   compute line btwn active corners based on corresponding active corner for midpoint??
    # for each interval
    #   compute g(x) based on bounds in dict
    #   compute line?? either here or above
    #   generate clause w this info

    # Generating boolean formula over (x,y) for unsafe region
    unsafe_cond = None  # initialize to None and build up with Ors
    for (traj1, traj2) in trajs:
        if unsafe_cond is None:
            unsafe_cond = traj1 * traj2 <= 0
        else:
            unsafe_cond = unsafe_cond | (traj1 * traj2 <= 0)
    if add_notches:
        for transition_point in set_of_transitions:
            # neg for left side, 0 for on edge, pos for on right side
            # inside polygon IF all neg or IF all pos
            shifted_vertex_pairs = [
                (v + transition_point, nextv + transition_point)
                for (v, nextv) in vertex_pairs
            ]  # [(v1, v2), (v2, v3), ..., (vn, v1)]
            # source: https://inginious.org/course/competitive-programming/geometry-pointinconvex
            # source: https://www.eecs.umich.edu/courses/eecs380/HANDOUTS/PROJ2/InsidePoly.html
            side_conds = [
                (y - v.y) * (nextv.x - v.x) - (x - v.x) * (nextv.y - v.y)
                for (v, nextv) in shifted_vertex_pairs
            ]  # [-1, 0, 4, 6]
            # construct both cases inside
            cond1 = true  # init to sympy.true since we're cascading Ands
            cond2 = true  # init to sympy.true since we're cascading Ands
            for side_cond in side_conds:
                cond1 = cond1 & (side_cond <= 0)
                cond2 = cond2 & (side_cond >= 0)
            unsafe_cond = unsafe_cond | cond1 | cond2

    # Translating to mathematica here for one-click plotting
    # boolean safe region condition for RegionPlot
    cond_mathematica: str = sympy_to_mathematica(unsafe_cond)
    # range over which to plot
    plotrange_mathematica: str = (
        f"{{x, {xbounds[0]}, {xbounds[1]}}}, {{y, {ybounds[0]}, {ybounds[1]}}}"
    )
    # trajectory equation to draw in Mathematica
    traj_eqn = Eq(trajectory, 0)
    # cut beginning 3 chars "Eq(" and last character ")"
    traj_mathematica: str = sympy_to_mathematica(str(traj_eqn)[3:-1]).replace(
        ", ", " == "
    )

    # TODO(nishant): add parameter or logic for offsetx
    offsetx = 1
    offsety = solve(traj_eqn.subs(x, offsetx))[0]
    mathematica_vertices = str([(v.x + offsetx, v.y + offsety) for v in verts])
    mathematica_vertices = "{" + sympy_to_mathematica(mathematica_vertices)[1:-1] + "}"
    mathematica_vertices = (
        "Polygon[" + mathematica_vertices.replace("(", "{").replace(")", "}") + "]"
    )

    TRAJ_COLOR = "Purple"
    POLY_COLOR = "Green"
    mathematica_header = (
        "\n======== Output for Plotting Safe Region in Mathematica ========\n"
    )
    mathematica_footer = (
        "\n=========================== Copy me! ==========================="
    )
    mathematica_output = f"""\nShow[
    RegionPlot[{cond_mathematica},  {plotrange_mathematica}, PlotPoints -> 60,  MaxRecursion -> 5],
    ContourPlot[{traj_mathematica},  {plotrange_mathematica}, ContourStyle->{{{TRAJ_COLOR}, Dashed}}],
    Graphics[ {{FaceForm[None], EdgeForm[{POLY_COLOR}], {mathematica_vertices} }} ],
    GridLines->Automatic, Ticks->Automatic\n]\n"""

    print(mathematica_header, mathematica_output, mathematica_footer)

    nelem = (xbounds[1] - xbounds[0]) * (ybounds[1] - ybounds[0]) / (resolution ** 2)
    count = 0
    for x0 in np.arange(xbounds[0], xbounds[1], resolution):
        for y0 in np.arange(ybounds[0], ybounds[1], resolution):
            count += 1
            # TODO(nishant): progress for plotting dots
            # print(f"{count/nelem*100}% \r")
            intruder = Point(x0, y0)
            is_safe = True
            for (traj1, traj2) in trajs:
                if (traj1.subs(x, intruder[0]).subs(y, intruder[1])) * (
                    traj2.subs(x, intruder[0]).subs(y, intruder[1])
                ) <= 0:
                    is_safe = False
                    break
            if not is_safe:
                for transition_point in set_of_transitions:
                    if is_safe and not encloses_method(
                        poly, transition_point, intruder
                    ):
                        is_safe = False
                        break
            if is_safe:
                ax.plot(x0, y0, "bo", alpha=alpha)
            else:
                ax.plot(x0, y0, "ro", alpha=alpha)

    ax.axis("equal")
    # p1 = plot_implicit(trajectory, line_color='k')
    # move_sympyplot_to_axes(p1, ax)

    ax.set_title(title)
    if savefig:
        plt.savefig(title)
    plt.show()
