import matplotlib.pyplot as plt
import numpy as np
import sympy

from sympy import *
from typing import Tuple, List, Dict, Set

init_printing(use_unicode=True)


def plot_polygon(poly: sympy.Polygon):
    # Draw a polygon by plotting vertices as points and edges as lines
    fig = plt.figure()
    ax = fig.gca()

    verts = poly.vertices

    for p in verts:
        ax.scatter(p.x, p.y, c="r")
        plt.text(p.x + 0.05, p.y + 0.05, f"({p.x},{p.y})")

    for (p, nextp) in zip(verts, verts[1:] + verts[:1]):
        x = np.linspace(float(p.x), float(nextp.x), 100, dtype=float)
        y = np.linspace(float(p.y), float(nextp.y), 100, dtype=float)
        ax.plot(x, y, c="b")

    ax.axis("equal")

    plt.show()


def compute_polygon_angles(poly: sympy.Polygon) -> Tuple[List, List]:
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


def find_transitions(trajectory, angles, x, y, domain=S.Complexes) -> Tuple[Dict, Set]:
    """
    TODO(nishant): add docstrings
    """
    # NOTE: trajectory is an *expression*, not equation
    transitions = {}
    # TODO(nishant): explain slope and variables for non f(x,y) functions
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


# TODO(nishant): this is broken, delete
def get_angle(
    theta, angles
):  # angles correspond to polygon edges, theta is the slope at a particular point
    theta = theta % (2 * pi)  # ensure theta is in the range [0, 2*pi)
    min_angle = min(angles)
    max_angle = max(angles)
    angle_pairs = list(zip(angles, angles[1:] + angles[:1]))
    if theta < min_angle or theta > max_angle:
        return (max_angle, min_angle)
    for (angle1, angle2) in angle_pairs:
        if theta == angle1:
            return angle1
        elif angle1 < theta and theta < angle2:
            return (angle1, angle2)
    assert False  # should never exit for loop without returning something


def outside_active_corners(
    poly: sympy.Polygon,
    trajectory,
    intruder: sympy.Point,
    angles_to_vertices: Dict,
    angle_range,
    x,
    y,
):
    assert (
        angles_to_vertices[angle_range[0]][1] == angles_to_vertices[angle_range[1]][0]
    )
    point1 = angles_to_vertices[angle_range[0]][1]
    point2 = poly.center - point1
    offset1 = point1 - poly.center
    offset2 = point2 - poly.center
    traj1 = trajectory.subs(x, x - offset1[0]).subs(y, y - offset1[1])
    traj2 = trajectory.subs(x, x - offset2[0]).subs(y, y - offset2[1])
    return (traj1.subs(x, intruder[0]).subs(y, intruder[1])) * (
        traj2.subs(x, intruder[0]).subs(y, intruder[1])
    ) > 0


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


# THIS IS BAD
# def safe(
#     poly: sympy.Polygon, trajectory, location: sympy.Point, intruder: sympy.Point, x, y
# ):
#     angles, vertex_pairs = compute_polygon_angles(poly)
#     angle_range = get_angle(atan2(*eval_slope(trajectory, location, x, y)), angles)
#     if type(angle_range) == tuple:  # use active corners
#         angles_to_vertices: Dict = dict(zip(angles, vertex_pairs))
#         return outside_active_corners(
#             poly, trajectory, intruder, angles_to_vertices, angle_range, x, y
#         )
#     else:  # check to see if inside polygon
#         return encloses_method(poly, location, intruder)


def move_sympyplot_to_axes(p, ax):
    backend = p.backend(p)
    backend.ax = ax
    # Fix for > sympy v1.5
    backend._process_series(backend.parent._series, ax, backend.parent)
    backend.ax.spines["right"].set_color("none")
    backend.ax.spines["bottom"].set_position("zero")
    backend.ax.spines["top"].set_color("none")
    plt.close(backend.fig)


# finding close parens
def find_parens(s: str) -> Dict[int, int]:
    # source: https://stackoverflow.com/a/29992065
    toret = {}
    pstack = []

    for i, c in enumerate(s):
        if c == "(":
            pstack.append(i)
        elif c == ")":
            if len(pstack) == 0:
                raise IndexError("No matching closing parens at: " + str(i))
            toret[pstack.pop()] = i

    if len(pstack) > 0:
        raise IndexError("No matching opening parens at: " + str(pstack.pop()))

    return toret


def swap_fn(s: str, swaps: List[Tuple[str, str]]) -> str:
    parens = find_parens(s)
    new = list(s)
    newstr = s
    for (tofind, toreplace) in swaps:
        while newstr.find(tofind) > -1:
            lendiff = len(toreplace) - len(tofind)
            first_idx = newstr.find(tofind)
            new = list(newstr)
            new[first_idx : first_idx + len(tofind)] = toreplace
            open_idx = newstr.find("(", first_idx)
            # account for different operator length
            new[open_idx + lendiff] = "["
            close_idx = parens[open_idx]
            # account for different operator length
            new[close_idx + lendiff] = "]"
            newstr = "".join(new)
    return newstr


def swap_piecewise(s: str):
    newstr = s
    while newstr.find("Piecewise(") > -1:
        parens = find_parens(newstr)
        val = newstr.find("Piecewise(") + len("Piecewise(") - 1
        newstr = newstr[:val] + "[{" + newstr[val + 1 :]
        newstr = newstr[: parens[val] + 1] + "}]" + newstr[parens[val] + 2 :]
        a = newstr[val + 2 : parens[val] + 1]  # input to Piecewise
        newa = []
        count = 0
        for (i, c) in enumerate(a):
            if c == "(":
                count += 1
                if count == 1:
                    newa.append("{")
                else:
                    newa.append("(")
            elif c == ")":
                count -= 1
                if count == 0:
                    newa.append("}")
                else:
                    newa.append(")")
            else:
                newa.append(c)

        a = "".join(newa)
        newstr = newstr[: val + 2] + a + newstr[parens[val] + 1 :]
    return newstr


# TODO(nishant): section for all the text transformations
def sympy_to_mathematica(cond) -> str:
    cond_string: str = str(cond)
    replacements = [
        ("**", "^"),
        ("|", "||"),
        ("&", "&&"),
        ("pi", "Pi"),
    ]
    for swap in replacements:
        assert (
            len(swap) == 2
        ), "There must be only two elements for specifying replacements"
        cond_string = cond_string.replace(swap[0], swap[1])

    functions = [
        ("sin", "Sin"),
        ("cos", "Cos"),
        ("tan", "Tan"),
        ("sqrt", "Sqrt"),
        # ("exp", "E"),  # maybe not the case
    ]

    cond_string = swap_fn(cond_string, functions)

    cond_string = swap_piecewise(cond_string)

    # TODO(nishant): return entire string to send to Mathematica
    #                with pred and RegionPlot[] call, domain, etc

    return cond_string


def plot_safe_grid_piecewise(
    x,
    y,
    poly: sympy.Polygon,
    trajectory,  # piecewise
    xbounds,
    ybounds,
    title,
    domain,
    resolution=0.25,
    alpha=1,
    savefig=False,
    add_notches=True,
    overlay_traj=False,
):
    """Generates clauses
    TODO(nishant): docstring
    """

    # TODO(nishant): figure out input spec
    # TODO(nishant): merge two ways of doing piecewise? prob not rn
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
            elif x not in subtraj.free_symbols:  # form x=f(y)
                _, subset_of_transitions = find_transitions(
                    -x + subtraj, angles, x, y, domain=subdomain
                )
            set_of_transitions.update(subset_of_transitions)

            # add piecewise boundary
            left_bound = Point(subdomain.inf, subtraj.subs(func_var, subdomain.inf))
            right_bound = Point(subdomain.sup, subtraj.subs(func_var, subdomain.sup))
            set_of_transitions.add(left_bound)
            set_of_transitions.add(right_bound)
    else:
        if y not in trajectory.free_symbols:
            _, subset_of_transitions = find_transitions(
                -y + trajectory, angles, x, y, domain=domain
            )
        elif x not in trajectory.free_symbols:
            _, subset_of_transitions = find_transitions(
                -x + trajectory, angles, x, y, domain=domain
            )

        set_of_transitions.update(subset_of_transitions)

        left_bound = Point(domain.inf, trajectory.subs(x, domain.inf))
        right_bound = Point(domain.sup, trajectory.subs(x, domain.sup))
        set_of_transitions.add(left_bound)
        set_of_transitions.add(right_bound)

    # TODO(nishant): get this working for f(y) also
    sorted_transitions: list = sorted(set_of_transitions, key=lambda point: point.x)
    x_transitions = [p.x for p in sorted_transitions]
    midpoints = np.convolve(x_transitions, [1, 1], "valid") / 2
    deriv_traj = diff(trajectory, x)
    dydx_midpoints = [deriv_traj.subs(x, val) for val in midpoints]
    midpoint_angles = [atan2(d, 1) for d in dydx_midpoints]
    # above needs to be modified to work for f(y) also!

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
    x_intervals = list(
        zip(x_transitions[:-1], x_transitions[1:])
    )  # same len as midpoints

    total_cond = None
    for i, x_interval in enumerate(x_intervals):
        x_left = x_interval[0]
        y_left = trajectory.subs(x, x_left)
        x_right = x_interval[1]
        y_right = trajectory.subs(x, x_right)
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
        x_cond = (x >= (x_left - w)) & (x <= (x_right + w))

        # construct g functions
        # trajectory over this interval, held constant outside of it
        if type(trajectory) == Piecewise:
            current_piece = trajectory.as_expr_set_pairs(Interval(x_left, x_right))
            current_fn = current_piece[0][0]
            g = y - Piecewise(
                (y_left, x < x_left), (current_fn, x <= x_right), (y_right, x > x_right)
            )
        else:
            g = y - Piecewise(
                (y_left, x < x_left), (trajectory, x <= x_right), (y_right, x > x_right)
            )
        #     TODO(nishant): rewrite for asymmetric
        g1 = g.subs(x, x - active_corner_offset.x).subs(y, y - active_corner_offset.y)
        g2 = g.subs(x, x + active_corner_offset.x).subs(y, y + active_corner_offset.y)

        full_cond = corner_cond & x_cond & (g1 * g2 <= 0)
        #     full_cond = (g1*g2 <= 0)
        #     print(full_cond)
        if total_cond is None:
            total_cond = full_cond
        else:
            total_cond = total_cond | full_cond
    #     print("")

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

    # Translating to mathematica here for one-click plotting
    # boolean safe region condition for RegionPlot
    cond_mathematica: str = sympy_to_mathematica(total_cond)
    # range over which to plot
    # TODO(nishant): fix bounds to be only one variable
    plotrange_mathematica: str = (
        f"{{x, {xbounds[0]}, {xbounds[1]}}}, {{y, {ybounds[0]}, {ybounds[1]}}}"
    )
    # trajectory equation to draw in Mathematica
    # TODO(nishant): prob not 0 now, maybe set equal to y?
    # traj_eqn = Eq(trajectory, 0)
    # cut beginning 3 chars "Eq(" and last character ")"
    # traj_mathematica: str = sympy_to_mathematica(str(traj_eqn)[3:-1]).replace(
    #     ", ", " == "
    # )
    traj_mathematica: str = sympy_to_mathematica(trajectory)

    # TODO(nishant): add parameter or logic for offsetx
    offsetx = 1
    # offsety = solve(traj_eqn.subs(x, offsetx))[0]
    offsety = trajectory.subs(x, offsetx)
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
        "\n=========================== Copy me! ===========================\n"
    )
    mathematica_output = f"""\nShow[
    RegionPlot[{cond_mathematica},  {plotrange_mathematica}, PlotPoints -> 60,  MaxRecursion -> 5],
    Plot[{traj_mathematica},  {plotrange_mathematica}, PlotStyle->{{{TRAJ_COLOR}, Dashed}}],
    Graphics[ {{FaceForm[None], EdgeForm[{POLY_COLOR}], {mathematica_vertices} }} ],
    GridLines->Automatic, Ticks->Automatic, AspectRatio->Automatic\n]\n"""

    print(mathematica_header, mathematica_output, mathematica_footer)

    if overlay_traj:
        p1 = plot(trajectory)

    fig = plt.figure(figsize=(15, 5))
    ax = fig.gca()

    if overlay_traj:
        move_sympyplot_to_axes(p1, ax)

    nelem = (xbounds[1] - xbounds[0]) * (ybounds[1] - ybounds[0]) / (resolution ** 2)
    count = 0
    for x0 in np.arange(xbounds[0], xbounds[1], resolution):
        for y0 in np.arange(ybounds[0], ybounds[1], resolution):
            count += 1
            # TODO(nishant): progress for plotting dots
            # print(f"{count/nelem*100}% \r")
            # intruder = Point(x0, y0)
            is_safe = (~total_cond).subs([(x, x0), (y, y0)])
            # is_safe = True
            # for (traj1, traj2) in trajs:
            #     if (traj1.subs(x, intruder[0]).subs(y, intruder[1])) * (
            #         traj2.subs(x, intruder[0]).subs(y, intruder[1])
            #     ) <= 0:
            #         is_safe = False
            #         break
            # if not is_safe:
            #     for transition_point in set_of_transitions:
            #         if is_safe and not encloses_method(
            #             poly, transition_point, intruder
            #         ):
            #             is_safe = False
            #             break
            # TODO(nishant): resolution scaling
            if is_safe:
                ax.plot(x0, y0, "bo", alpha=alpha, markersize=resolution * 4)
            else:
                ax.plot(x0, y0, "ro", alpha=alpha, markersize=resolution * 4)

    ax.axis("equal")
    # p1 = plot_implicit(trajectory, line_color='k')
    # move_sympyplot_to_axes(p1, ax)

    ax.set_title(title)
    if savefig:
        plt.savefig(title)
    plt.show()


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


# traj -> [x(t); y(t)] -> at some T, what is the angle of the tangent to trajectory
# x,y points, may or may not be on the trajectory,

# at any given (x,y) on traj, there's either a notch OR some pair is active,
# and bounds everything else -> suffices to check for notches AND to check
# that the point is outside *all* active-corner pairs

# TODO(nishant): write better in spec
# given obstacle at (x_O, y_O)
# algorithm is check for polygon inclusion at transition points AND check point
# is outside all active-corner-pairs

# assume symmetric polygons and write something to identify the active corner pairs
# try with a hexagon, square, diamond, rectangle


# Implemented the below two functions when I was confused about how to
# identify active corners in between transition points - but we can just
# test safety using all active-corner pairs so this stuff isn't required.


# def find_common_corner(angles_to_vertices: dict, angle, direction):
#     # I don't think this is actually useful but keeping it around
#     assert abs(direction) == 1
#
#     sorted_angles = sorted(angles_to_vertices.keys())
#
#     angle_index = sorted_angles.index(angle)
#     next_greatest_angle_index = (angle_index + direction) % len(
#         angles
#     )  # wrap around 2pi
#     next_greatest_angle = sorted_angles[next_greatest_angle_index]
#
#     common_corner: set = set(angles_to_vertices[angle]).intersection(
#         set(angles_to_vertices[next_greatest_angle])
#     )
#     assert len(common_corner) == 1
#     common_corner: Point = common_corner.pop()
#     return common_corner
#
#
# def direction_of_traj_angle(trajectory, transition_points, angle, epsilon):
#     # I don't think this is actually useful but keeping it around
#     assert epsilon != 0, "Epsilon cannot be 0"
#
#     p = transition_points[angle][0]  # TODO: fix hard coding
#     m = tan(angle)
#     y_on_tangent = p.y + epsilon * m
#     y_solns = solve(Eq(trajectory.subs(x, p.x + epsilon), 0))  # two solutions for y
#
#     # find closest to y_on_tangent
#     closest_y_soln = min(y_solns, key=lambda v: abs(v - y_on_tangent))
#
#     if epsilon > 0:
#         if closest_y_soln > y_on_tangent:
#             # theta increases: pick corners accordingly
#             direction = +1
#         else:
#             # theta decreases: pick corners accordingly
#             direction = -1
#     elif epsilon < 0:
#         if closest_y_soln > y_on_tangent:
#             # theta decreases: pick corners accordingly
#             direction = -1
#         else:
#             # theta increases: pick corners accordingly
#             direction = +1
#     return direction
