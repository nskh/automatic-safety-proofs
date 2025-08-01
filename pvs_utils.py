import numpy as np
import matplotlib.pyplot as plt
from sympy import expand, Line, Polygon, symbols, Point, Function, diff, atan2, pi, oo
from sympy.functions.elementary.piecewise import Piecewise


# =============================================================================
# PLOTTING AND VISUALIZATION UTILITIES
# =============================================================================


def plot_polygon(poly: Polygon, labels=[]):
    """Plot a polygon with vertices and edges."""
    fig = plt.figure()
    ax = fig.gca()

    verts = poly.vertices

    for i, p in enumerate(verts):
        ax.scatter(p.x, p.y, c="r")
        if labels:
            plt.text(p.x + 0.05, p.y + 0.05, f"{labels[i]}: ({p.x},{p.y})")
        else:
            plt.text(p.x + 0.05, p.y + 0.05, f"({p.x},{p.y})")

    for p, nextp in zip(verts, verts[1:] + verts[:1]):
        x = np.linspace(float(p.x), float(nextp.x), 100, dtype=float)
        y = np.linspace(float(p.y), float(nextp.y), 100, dtype=float)
        ax.plot(x, y, c="b")

    ax.axis("equal")
    plt.show()


# =============================================================================
# PREMISE AND CONDITION GENERATION
# =============================================================================


def generate_premise(lines: dict, traj_expr):
    """Generate a premise from lines and trajectory expression."""
    lines_rel = make_relative(lines, traj_expr)
    premise = True
    for line_rel in lines_rel.values():
        premise &= line_rel
    return premise


def generate_explicit(
    corner_pairs: list,
    traj,
    verts,
    x=symbols("x"),
    y=symbols("y"),
    xo=symbols("xo"),
    yo=symbols("yo"),
):
    """Generate explicit conditions for corner pairs."""
    explicit = False
    for c1, c2 in corner_pairs:
        corner_traj_1 = expand(traj.subs(y, yo - verts[c1].y).subs(x, xo - verts[c1].x))
        corner_traj_2 = expand(traj.subs(y, yo - verts[c2].y).subs(x, xo - verts[c2].x))
        pair_clause = corner_traj_1 * corner_traj_2 <= 0
        explicit |= pair_clause
    return explicit


def generate_explicit_disjunction(
    corner_pairs: list,
    traj,
    verts,
    x=symbols("x"),
    y=symbols("y"),
    xo=symbols("xo"),
    yo=symbols("yo"),
):
    """Generate explicit disjunction conditions for corner pairs."""
    explicit = False
    for c1, c2 in corner_pairs:
        corner_traj_1 = expand(traj.subs(y, yo - verts[c1].y).subs(x, xo - verts[c1].x))
        corner_traj_2 = expand(traj.subs(y, yo - verts[c2].y).subs(x, xo - verts[c2].x))
        pair_clause = ((corner_traj_1 <= 0) & (corner_traj_2 >= 0)) | (
            (corner_traj_1 >= 0) & (corner_traj_2 <= 0)
        )
        explicit |= pair_clause
    return explicit


# =============================================================================
# PVS FORMATTING AND CONVERSION
# =============================================================================


def sympy_to_pvs(clause: str):
    """Convert sympy expressions to PVS format."""
    if type(clause) is not str:
        clause = str(clause)
    return (
        clause.replace("x_O", "xo")
        .replace("y_O", "yo")
        .replace("|", "OR\n   ")
        .replace("&", "AND\n    ")
        .replace("**", "^")
    )


def construct_lemma_old(premise, active_corner_condition, lemma_name):
    """Legacy lemma construction function."""
    return f"""{lemma_name}(xo, yo, alpha: real) : bool =
    (EXISTS (x : real) :
    ({sympy_to_pvs(str(premise))})
    IMPLIES
    {sympy_to_pvs(str(active_corner_condition))}
    """


def full_lemma(
    poly,
    vert_names,
    corner_pairs,
    traj,
    traj_expr=None,
    lemma_name="Soundness",
    disjunction=False,
):
    """Generate a full lemma with polygon and trajectory information."""
    if type(vert_names) is str:
        vert_names = vert_names.split()
    assert type(vert_names) is list
    verts, lines = verts_and_lines(vert_names, poly)

    premise = generate_premise(lines, traj_expr)
    if disjunction:
        explicit = generate_explicit_disjunction(corner_pairs, traj, verts)
    else:
        explicit = generate_explicit(corner_pairs, traj, verts)
    return construct_lemma_old(premise, explicit, lemma_name)


# =============================================================================
# GEOMETRIC UTILITIES
# =============================================================================


def build_line(x, y, coefs):
    """Build a line equation from coefficients."""
    return x * coefs[0] + y * coefs[1] + coefs[2]


def make_relative(
    expr,
    traj_expr,
    x=symbols("x"),
    y=symbols("y"),
    xo=symbols("xo"),
    yo=symbols("yo"),
    alpha=symbols("alpha"),
):
    """Make expressions relative to trajectory."""
    if type(expr) == dict:
        return make_relative_dict(expr, traj_expr)
    if traj_expr is None:
        return expr.subs(x, xo - x).subs(y, yo - alpha * x)
    else:
        return expr.subs(x, xo - x).subs(y, yo - traj_expr)


def make_relative_dict(d, te):
    """Make dictionary of expressions relative to trajectory."""
    return {k: make_relative(v, te) for k, v in d.items()}


def verts_and_lines(
    vert_names: list[str], poly: Polygon, x=symbols("x"), y=symbols("y")
):
    """Extract vertices and lines from polygon."""
    assert len(vert_names) == len(poly.vertices)
    verts: dict = dict(zip(vert_names, poly.vertices))
    vert_names = sorted(vert_names)
    vert_pairs: list = list(zip(vert_names, vert_names[1:] + vert_names[0:1]))
    lines: dict = {}
    for i, (vert1, vert2) in enumerate(vert_pairs, 1):
        line = build_line(x, y, Line(verts[vert1], verts[vert2]).coefficients)
        # sub in (0,0) to find inequality direction
        if (line >= 0).subs([(x, 0), (y, 0)]):
            lines[i] = line >= 0
        else:
            lines[i] = line <= 0
    return verts, lines


# =============================================================================
# LEMMA CONSTRUCTION
# =============================================================================


def construct_lemma(
    verts,
    active_exists_premise,
    active_corner_condition,
    lemma_name,
    deriv_clause1=" >= 0)",
    deriv_clause2=None,
):
    """
    Produces a lemma string of the form:

      lemma_name: LEMMA
        FORALL(f:[real-> real],xo,yo:real):
        deriv_clause1 AND
        deriv_clause2 AND
        deriv_clause3 AND
        (EXISTS (x : real) : (active_exists_premise)
           AND exists_upper >= x AND exists_lower <= x)
        IMPLIES
        active_corner_condition
    """
    # Build the exists clause from the main exists-premise and the two bounds.
    max_offset = max([v.x for v in verts.values()])
    exists_upper = f"xo + {max_offset}"

    min_offset = min([v.x for v in verts.values()])
    exists_lower = f"xo + {min_offset}"
    exists_clause = f"""(EXISTS (x : real) :
    ({sympy_to_pvs(str(active_exists_premise))}) AND
    {exists_upper} >= x AND {exists_lower} <= x)"""

    differentiable_statement = "derivable?[real](f)"
    if deriv_clause2:
        full_preamble = f"""{differentiable_statement} AND
    (FORALL(x:real): deriv[real](f)(x) {deriv_clause1}) AND
    (FORALL(x:real): deriv[real](f)(x) {deriv_clause2}) AND
    {exists_clause}"""
    else:
        full_preamble = f"""{differentiable_statement} AND
    (FORALL(x:real): deriv[real](f)(x) {deriv_clause1}) AND
    {exists_clause}"""

    lemma_str = f"""{lemma_name}: LEMMA
    FORALL(f:[real-> real],xo,yo:real):
    {full_preamble}
    IMPLIES
    {sympy_to_pvs(str(active_corner_condition))}
"""
    return lemma_str


# =============================================================================
# PROOF NODE AND SCRIPT CLASSES
# =============================================================================


def unbounded_one_side_proof_script(
    case_label,
    deriv_lemma,
    max_right,
    min_left,
    domain_definition,
    max_right_clipped,
    min_left_clipped,
):
    """For domains unbounded to the left/right, we need to handle two cases:
    - Case 1: Unclipped trajectory
    - Case 2: Clipped trajectory

    Inputs:
    - case_label: The case label for the CASE node, should be TRUE for unclipped trajectory and FALSE for clipped trajectory
    - deriv_lemma: The name of the derivative lemma to use: mvt_gen_ge_lo_2, etc
    - max_right: The maximum right value inside polygon, e.g. xo+half_width
    - min_left: The minimum left value inside polygon, e.g. xo-half_width
    - domain_definition: The definition of the domain, e.g. "left_open" or "right_open"
    - max_right_clipped: The maximum right value after clipping, e.g. xo+half_width but maybe 0
    - min_left_clipped: The minimum left value after clipping, e.g. xo-half_width
    """
    return f"""%|- (THEN (SKEEP*) (SKOLETIN*) (FLATTEN) (SKEEP)
    %|-  (SPREAD (CASE "{case_label}")
    %|-   ((THEN (LEMMA "{deriv_lemma}")
    %|-     (SPREAD (INST -1 "f" "0" "0" "{max_right}" "x")
    %|-      ((SPREAD (SPLIT -1)
    %|-        ((THEN (ASSERT) (LEMMA "{deriv_lemma}")
    %|-          (SPREAD (INST -1 "f" "0" "0" "x" "{min_left}")
    %|-           ((ASSERT) (THEN (EXPAND "{domain_definition}") (ASSERT))
    %|-            (THEN (EXPAND "{domain_definition}") (ASSERT)))))
    %|-         (PROPAX) (PROPAX) (PROPAX)))
    %|-       (THEN (EXPAND "{domain_definition}") (ASSERT))
    %|-       (THEN (EXPAND "{domain_definition}") (ASSERT)))))
    %|-    (THEN (LEMMA "{deriv_lemma}")
    %|-     (SPREAD (INST -1 "f" "0" "0" "{max_right_clipped}" "x")
    %|-      ((SPREAD (SPLIT -1)
    %|-        ((THEN (ASSERT) (LEMMA "{deriv_lemma}")
    %|-          (SPREAD (INST -1 "f" "0" "0" "x" "{min_left_clipped}")
    %|-           ((THEN (EXPAND "f") (EXPAND "{domain_definition}") (SPREAD (SPLIT -1) ((ASSERT) (PROPAX))))
    %|-            (THEN (EXPAND "f") (EXPAND "{domain_definition}") (ASSERT))
    %|-            (THEN (EXPAND "f") (EXPAND "{domain_definition}") (ASSERT)))))
    %|-         (ASSERT) (PROPAX) (PROPAX)))
    %|-       (THEN (EXPAND "f") (EXPAND "{domain_definition}") (ASSERT))
    %|-       (THEN (EXPAND "f") (EXPAND "{domain_definition}") (ASSERT))))))))
    """


class ProofNode:
    """Represents a node in a PVS proof tree."""

    def __init__(self, command, args=None, kwargs=None, children=None):
        self.command = command
        self.args = args if args else []
        self.kwargs = kwargs if kwargs else {}
        self.children = children if children is not None else []

    def __repr__(self):
        args_str = f" {self.args}" if self.args else ""
        kwargs_str = f" {self.kwargs}" if self.kwargs else ""
        children_str = f" children={len(self.children)}" if self.children else ""
        return f"ProofNode({self.command}{args_str}{kwargs_str}{children_str})"

    def add_child(self, node):
        """Add a child node."""
        self.children.append(node)
        return node

    def add_children(self, nodes):
        """Add multiple child nodes."""
        self.children.extend(nodes)
        return nodes

    def generate(self, indent="%|- ", depth=0, n_spaces=4):
        """Generate PVS proof script from the node tree."""
        base_indent = indent + " " * depth * n_spaces
        # Compose argument string
        arg_str = " ".join(str(arg) for arg in self.args) if self.args else ""
        kwarg_str = (
            " ".join(f":{k} {v}" for k, v in self.kwargs.items()) if self.kwargs else ""
        )
        all_args = " ".join(filter(None, [arg_str, kwarg_str]))
        # Special handling for wrapper nodes
        if self.command == "()":
            child_proofs = [c.generate(indent, depth) for c in self.children]
            return f"({' '.join(child_proofs)})"
        # Base case - no children
        if not self.children:
            return f"({self.command}{' ' + all_args if all_args else ''})"
        # Multiple children or nested structures
        result = f"({self.command}{' ' + all_args if all_args else ''}"
        child_proofs = [c.generate(indent, depth + 1) for c in self.children]
        child_str = f"\n{base_indent} ".join(child_proofs)
        return f"{result}\n{base_indent} {child_str})"


class ProofScript:
    """Represents a complete PVS proof script."""

    def __init__(self, lemma_name):
        self.lemma_name = lemma_name
        self.root = None

    def generate(self):
        """Generate the complete proof script."""
        try:
            lines = []
            lines.append(f"%|- {self.lemma_name} : PROOF")
            lines.append(f"%|- {self.root.generate()}")
            lines.append(f"%|- QED {self.lemma_name}")
            return "\n".join(lines)
        except Exception as e:
            print(self.lemma_name)
            print(self.root)
            print(e)


# =============================================================================
# AUTOMATIC PROOF SCRIPT GENERATION
# =============================================================================


def generate_unbounded_proof_calls(
    trajectory, poly, domain, x=symbols("x"), y=symbols("y")
):
    """
    Automatically generate calls to unbounded_one_side_proof_script based on
    trajectory, polygon, and domain analysis.

    Args:
        trajectory: Sympy expression or Piecewise for trajectory
        poly: Sympy Polygon object
        domain: Sympy Interval domain
        x, y: Sympy symbols for variables

    Returns:
        list: List of dictionaries containing parameters for unbounded_one_side_proof_script calls
    """
    # Import required functions
    try:
        from safe_region_utils import (
            compute_polygon_angles,
            find_transitions,
            compute_corners_to_angles,
        )
    except ImportError:
        # Fallback if direct import fails
        import sys
        import os

        sys.path.append(os.path.dirname(os.path.abspath(__file__)))
        from safe_region_utils import (
            compute_polygon_angles,
            find_transitions,
            compute_corners_to_angles,
        )

    import numpy as np

    # Trajectory is assumed to be function of x only
    func_var = x

    # Compute polygon width using same method as construct_lemma
    max_offset = max([v.x for v in poly.vertices])
    min_offset = min([v.x for v in poly.vertices])

    # Calculate half-width (distance from centroid)
    half_width = (max_offset - min_offset) / 2

    # Get polygon angles and find transitions
    angles, vertex_pairs = compute_polygon_angles(poly)
    set_of_transitions = set()

    if isinstance(trajectory, Piecewise):
        # Handle piecewise trajectories
        for subtraj, subcond in trajectory.as_expr_set_pairs():
            subdomain = subcond.intersect(domain)

            # Trajectory is y = f(x)
            _, subset_of_transitions = find_transitions(
                -y + subtraj, angles, x, y, domain=subdomain
            )
            left_bound = Point(subdomain.inf, subtraj.subs(func_var, subdomain.inf))
            right_bound = Point(subdomain.sup, subtraj.subs(func_var, subdomain.sup))

            set_of_transitions.update(subset_of_transitions)

            if left_bound.x.is_finite and left_bound.y.is_finite:
                set_of_transitions.add(left_bound)
            if right_bound.x.is_finite and right_bound.y.is_finite:
                set_of_transitions.add(right_bound)
    else:
        # Handle single trajectory (y = f(x))
        _, subset_of_transitions = find_transitions(
            -y + trajectory, angles, x, y, domain=domain
        )
        # Add left and right boundaries to check for notch there too
        left_bound = Point(domain.inf, trajectory.subs(func_var, domain.inf))
        right_bound = Point(domain.sup, trajectory.subs(func_var, domain.sup))

        set_of_transitions.update(subset_of_transitions)
        if left_bound.x.is_finite and left_bound.y.is_finite:
            set_of_transitions.add(left_bound)
        if right_bound.x.is_finite and right_bound.y.is_finite:
            set_of_transitions.add(right_bound)

    # Sort transitions and compute midpoints
    sorted_transitions = sorted(set_of_transitions, key=lambda point: point.x)
    func_var_transitions = [p.x for p in sorted_transitions]
    midpoints = np.convolve(func_var_transitions, [1, 1], "valid") / 2

    # Compute derivative and angles at midpoints
    if isinstance(trajectory, Piecewise):
        # For piecewise, we need to handle each piece separately
        deriv_traj = trajectory.diff(x)
    else:
        deriv_traj = diff(trajectory, x)

    dydx_midpoints = []
    for val in midpoints:
        try:
            deriv_val = deriv_traj.subs(func_var, val)
            dydx_midpoints.append(deriv_val)
        except:
            # Handle cases where derivative evaluation fails
            dydx_midpoints.append(0)  # Default to 0 if evaluation fails
    midpoint_angles = [atan2(d, 1) for d in dydx_midpoints]

    # Map angles to active corners (using existing logic from safe_region_utils)
    corners_to_angles = compute_corners_to_angles(poly)
    active_corners = {}
    for midpoint_angle in midpoint_angles:
        for k, v in corners_to_angles.items():
            if midpoint_angle % (2 * pi) in v:
                active_corners[midpoint_angle] = k

    # Generate intervals between transition points and domain boundaries
    # Include domain boundaries in the transition points
    domain_min = domain.inf
    domain_max = domain.sup

    # Add domain boundaries to transition points
    all_transition_points = func_var_transitions.copy()
    if domain_min.is_finite:
        all_transition_points.insert(0, domain_min)
    elif domain_min == -oo:
        all_transition_points.insert(0, -oo)
    if domain_max.is_finite:
        all_transition_points.append(domain_max)
    elif domain_max == oo:
        all_transition_points.append(oo)

    # Create intervals between all transition points
    var_intervals = list(zip(all_transition_points[:-1], all_transition_points[1:]))

    proof_calls = []

    for i, var_interval in enumerate(var_intervals):
        interval_start, interval_end = var_interval

        # Skip if this interval is outside the domain
        if interval_end <= domain_min or interval_start >= domain_max:
            continue

        # Determine if this interval is unbounded on left or right
        left_unbounded = not interval_start.is_finite
        right_unbounded = not interval_end.is_finite

        if not (left_unbounded or right_unbounded):
            # For bounded intervals, we'll treat them as left-unbounded for now
            left_unbounded = True
            right_unbounded = False

        # Find the active corner for this interval
        # Use the midpoint of this interval
        interval_midpoint = (interval_start + interval_end) / 2
        deriv_at_midpoint = deriv_traj.subs(func_var, interval_midpoint)

        # Map to active corner using existing logic
        # Find which midpoint_angles is closest to our interval_midpoint
        closest_midpoint_idx = min(
            range(len(midpoints)), key=lambda j: abs(midpoints[j] - interval_midpoint)
        )
        active_corner_offset = active_corners[midpoint_angles[closest_midpoint_idx]]

        # Generate case label for unbounded domains
        if left_unbounded:  # Interval unbounded on left
            case_label = f"xo + {half_width} >= {interval_end}"
        else:  # right_unbounded: Interval unbounded on right
            case_label = f"xo - {half_width} <= {interval_start}"

        # Determine derivative lemma based on derivative sign at midpoint
        deriv_sign = "ge" if deriv_at_midpoint >= 0 else "le"
        domain_type = "lo" if left_unbounded else "ro"
        deriv_lemma = f"mvt_gen_{deriv_sign}_{domain_type}_2"

        # Polygon bounds
        max_right = f"xo + {half_width}"
        min_left = f"xo - {half_width}"

        # Domain definition
        domain_definition = "left_open" if left_unbounded else "right_open"

        # Clipped bounds
        if left_unbounded:
            max_right_clipped = f"min(xo + {half_width}, {interval_end})"
            min_left_clipped = f"max(xo - {half_width}, {interval_start})"
        else:  # right_unbounded
            max_right_clipped = f"min(xo + {half_width}, {interval_end})"
            min_left_clipped = f"max(xo - {half_width}, {interval_start})"

        proof_call = {
            "case_label": case_label,
            "deriv_lemma": deriv_lemma,
            "max_right": max_right,
            "min_left": min_left,
            "domain_definition": domain_definition,
            "max_right_clipped": max_right_clipped,
            "min_left_clipped": min_left_clipped,
            "interval": var_interval,
            "active_corner": active_corner_offset,
        }

        proof_calls.append(proof_call)

    return proof_calls


def generate_proof_scripts_from_calls(proof_calls):
    """
    Generate actual proof script strings from proof call parameters.

    Args:
        proof_calls: List of dictionaries from generate_unbounded_proof_calls

    Returns:
        list: List of proof script strings
    """
    proof_scripts = []

    for call_params in proof_calls:
        script = unbounded_one_side_proof_script(
            case_label=call_params["case_label"],
            deriv_lemma=call_params["deriv_lemma"],
            max_right=call_params["max_right"],
            min_left=call_params["min_left"],
            domain_definition=call_params["domain_definition"],
            max_right_clipped=call_params["max_right_clipped"],
            min_left_clipped=call_params["min_left_clipped"],
        )
        proof_scripts.append(script)

    return proof_scripts
