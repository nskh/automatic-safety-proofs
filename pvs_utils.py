import numpy as np
import matplotlib.pyplot as plt
from sympy import expand, Line, Polygon, symbols, Point, Function, diff, atan2, pi, oo
from sympy.functions.elementary.piecewise import Piecewise
import string
from safe_region_utils import (
    compute_polygon_angles,
    find_transitions,
    compute_corners_to_angles,
)


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


def clean_trajectory(trajectory: str):
    """Clean trajectory COND so it's one line."""
    return " ".join([term.strip() for term in trajectory.split("\n")])


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
    domain_start=None,
    domain_end=None,
    notch=None,
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
    # closed interval: [domain_start, domain_end]
    if domain_start is not None and domain_end is not None:
        function_statement = f"""LET f = LAMBDA(x:real):
        COND
            x < {domain_start} -> g({domain_start}),
            x >= {domain_start} AND x <= {domain_end} -> g(x),
            ELSE -> g({domain_end})
        ENDCOND"""

        deriv_domain = f"(ci({domain_start}, {domain_end}))"
        exists_domain = f"x >= {domain_start} AND x <= {domain_end}"
    # right open interval: [domain_start, oo)
    elif domain_start is not None:
        function_statement = f"""LET f = LAMBDA(x:real):
        COND
            x < {domain_start} -> g({domain_start}),
            ELSE -> g(x)
        ENDCOND"""

        deriv_domain = f"(right_open({domain_start}))"
        exists_domain = f"x >= {domain_start}"
    # left open interval: (-oo, domain_end)
    elif domain_end is not None:
        function_statement = f"""LET f = LAMBDA(x:real):
        COND
            x <= {domain_end} -> g(x),
            ELSE -> g({domain_end})
        ENDCOND"""
        deriv_domain = f"(left_open({domain_end}))"
        exists_domain = f"x <= {domain_end}"
    else:
        raise ValueError("No domain specified")

    clipped_trajectory = function_statement[function_statement.find("LAMBDA") :]

    # Build the exists clause from the main exists-premise and the two bounds.
    max_offset = max([v.x for v in verts.values()])
    exists_upper = f"xo + {max_offset}"

    min_offset = min([v.x for v in verts.values()])
    exists_lower = f"xo + {min_offset}"
    exists_clause = f"""(EXISTS (x : real) :
    ({sympy_to_pvs(str(active_exists_premise))}) AND
    {exists_domain} AND
    {exists_upper} >= x AND {exists_lower} <= x)"""

    if deriv_clause2:
        deriv_list = [
            f"derivable?[{deriv_domain}](f)",
            f"(FORALL(x:{deriv_domain}): deriv[{deriv_domain}](f)(x) {deriv_clause1})",
            f"(FORALL(x:{deriv_domain}): deriv[{deriv_domain}](f)(x) {deriv_clause2})",
        ]
        deriv_statement = " AND ".join(deriv_list) + " AND"
        full_preamble = f"""    {function_statement}
    IN 
    {deriv_statement}
    {exists_clause}"""
    else:
        deriv_list = [
            f"derivable?[{deriv_domain}](f)",
            f"(FORALL(x:{deriv_domain}): deriv[{deriv_domain}](f)(x) {deriv_clause1})",
        ]
        deriv_statement = " AND ".join(deriv_list) + " AND"
        full_preamble = f"""    {function_statement}
    IN 
    {deriv_statement}
    {exists_clause}"""

    lemma_str = f"""{lemma_name}: LEMMA
    FORALL(xo,yo:real, g:[real -> real]):
    {full_preamble}
    IMPLIES
    {sympy_to_pvs(str(active_corner_condition))} OR
    {sympy_to_pvs(str(notch))}
"""
    return lemma_str, deriv_list, clipped_trajectory


def generate_unifying_lemma(
    verts,
    active_exists_premise,
    active_corner_condition,
    deriv_statements=[],
    notches=[],
    trajectories=[],
):
    max_offset = max([v.x for v in verts.values()])
    exists_upper = f"xo + {max_offset}"

    min_offset = min([v.x for v in verts.values()])
    exists_lower = f"xo + {min_offset}"
    exists_clause = f"""(({sympy_to_pvs(str(active_exists_premise))}) AND
    {exists_upper} >= x AND {exists_lower} <= x)"""

    notch_string = (
        " OR\n".join([f"    ({sympy_to_pvs(str(notch))})" for notch in notches])
        .replace("(f)", "(g)")
        .replace("f(", "g(")
    )

    deriv_premise = " AND\n    ".join(deriv_statements)
    # Stitch all the trajectories together and also represent the active corner condition with each of the clipped functions `f{i}`
    trajectory_statement = "LET\n"
    conditions = ""
    for i, trajectory in enumerate(trajectories):
        # do not add trailing "," for last case
        if i < len(trajectories) - 1:
            trajectory_statement += f"        f{i} = {trajectory},\n"
        else:
            trajectory_statement += f"        f{i} = {trajectory}\n"
        explicit_condition = sympy_to_pvs(str(active_corner_condition)).replace(
            "f(", f"f{i}("
        )
        conditions += explicit_condition + " OR\n"

    print(trajectory_statement)

    full_preamble = f"""{trajectory_statement}
    IN
    {deriv_premise} AND
    {exists_clause}"""

    lemma_text = f"""full_domain_soundness_lemma: LEMMA
    FORALL(x,xo,yo:real, g:[real -> real]):
    {full_preamble.replace("(f)", "(g)").replace("f(", "g(")}
    IMPLIES
    {conditions}{notch_string}
"""
    return lemma_text


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
    domain_bound,
    deriv_bound,
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
    - domain_bound: The domain bound for things like left_open, right_open, ci, e.g. "0"
    - deriv_bound: The derivative bound, e.g. "0"
    """

    # trim the (bound) off of domain definitions
    domain_definition = domain_definition.split("(")[0]

    return f"""%|- (THEN (SKEEP*) (SKOLETIN*) (FLATTEN) (SKEEP)
%|-  (SPREAD (CASE "{case_label}")
%|-   ((THEN (LEMMA "{deriv_lemma}")
%|-     (SPREAD (INST -1 "f" "{domain_bound}" "{deriv_bound}" "{max_right}" "x")
%|-      ((SPREAD (SPLIT -1)
%|-        ((THEN (ASSERT) (LEMMA "{deriv_lemma}")
%|-          (SPREAD (INST -1 "f" "{domain_bound}" "{deriv_bound}" "x" "{min_left}")
%|-           ((ASSERT) (THEN (EXPAND "{domain_definition}") (ASSERT))
%|-            (THEN (EXPAND "{domain_definition}") (ASSERT)))))
%|-         (PROPAX) (PROPAX) (PROPAX)))
%|-       (THEN (EXPAND "{domain_definition}") (ASSERT))
%|-       (THEN (EXPAND "{domain_definition}") (ASSERT)))))
%|-    (THEN (LEMMA "{deriv_lemma}")
%|-     (SPREAD (INST -1 "f" "{domain_bound}" "{deriv_bound}" "{max_right_clipped}" "x")
%|-      ((SPREAD (SPLIT -1)
%|-        ((THEN (ASSERT) (LEMMA "{deriv_lemma}")
%|-          (SPREAD (INST -1 "f" "{domain_bound}" "{deriv_bound}" "x" "{min_left_clipped}")
%|-           ((THEN (EXPAND "f") (EXPAND "{domain_definition}") (SPREAD (SPLIT -1) ((ASSERT) (PROPAX))))
%|-            (THEN (EXPAND "f") (EXPAND "{domain_definition}") (ASSERT))
%|-            (THEN (EXPAND "f") (EXPAND "{domain_definition}") (ASSERT)))))
%|-         (THEN (ASSERT) (PROPAX)) (THEN (ASSERT) (PROPAX)) (THEN (ASSERT) (PROPAX))))
%|-       (THEN (EXPAND "f") (EXPAND "{domain_definition}") (ASSERT))
%|-       (THEN (EXPAND "f") (EXPAND "{domain_definition}") (ASSERT))))))))
"""


def bounded_proof_script(
    case_label1,
    case_label2,
    case_label3,
    case_label4,
    deriv_lemma,
    max_right,
    min_left,
    max_right_clipped,
    min_left_clipped,
    domain_start,
    domain_end,
    deriv_bound,
):
    """For bounded domains, we need to handle four cases:
    - Case 1: Main case, in between lower and upper bounds of domain
    - Case 2: Outside domain, below lower bound
    - Case 3: Outside domain, above upper bound
    - Case 4: Odd placeholder case, below left but above right, we just assert out of this

    Inputs:
    - case_label1
    - case_label2
    - case_label3
    - case_label4
    - deriv_lemma: The name of the derivative lemma to use: mvt_gen_ge_ci, mvt_gen_le_ci
    - max_right: The maximum right value inside polygon, e.g. xo+half_width
    - min_left: The minimum left value inside polygon, e.g. xo-half_width
    - max_right_clipped: The maximum right value after clipping, e.g. xo+half_width but maybe 0
    - min_left_clipped: The minimum left value after clipping, e.g. xo-half_width

    Note that domain_definition is not passed in as we always use "ci".
    """
    return f"""%|- (THEN (SKEEP*) (SKOLETIN*) (FLATTEN) (SKEEP)
%|-  (SPREAD (CASE "{case_label1}")
%|-   ((THEN (FLATTEN) (LEMMA "{deriv_lemma}")
%|-     (SPREAD (INST -1 "f" "{domain_start}" "{domain_end}" "{deriv_bound}" "{max_right}" "x")
%|-      ((SPREAD (SPLIT -1)
%|-        ((THEN (ASSERT) (LEMMA "{deriv_lemma}")
%|-          (SPREAD (INST -1 "f" "{domain_start}" "{domain_end}" "{deriv_bound}" "x" "{min_left}")
%|-           ((ASSERT) (THEN (EXPAND "ci") (PROPAX))
%|-            (THEN (ASSERT) (EXPAND "ci") (PROPAX)))))
%|-         (ASSERT) (PROPAX) (PROPAX) (PROPAX)))
%|-       (THEN (EXPAND "ci") (ASSERT)) (THEN (EXPAND "ci") (ASSERT)))))
%|-    (SPREAD (CASE "{case_label2}")
%|-     ((THEN (FLATTEN) (LEMMA "{deriv_lemma}")
%|-       (SPREAD (INST -1 "f" "{domain_start}" "{domain_end}" "{deriv_bound}" "{max_right}" "x")
%|-        ((SPREAD (SPLIT -1)
%|-          ((THEN (ASSERT) (LEMMA "{deriv_lemma}")
%|-            (SPREAD (INST -1 "f" "{domain_start}" "{domain_end}" "{deriv_bound}" "x" "{min_left_clipped}")
%|-             ((THEN (EXPAND "f") (ASSERT)) (THEN (EXPAND "ci") (PROPAX))
%|-              (THEN (ASSERT) (EXPAND "ci") (PROPAX)))))
%|-           (ASSERT) (PROPAX) (PROPAX) (PROPAX)))
%|-         (THEN (EXPAND "ci") (ASSERT)) (THEN (EXPAND "ci") (ASSERT)))))
%|-      (THEN (ASSERT)
%|-       (SPREAD (CASE "{case_label3}")
%|-        ((THEN (FLATTEN) (LEMMA "{deriv_lemma}")
%|-          (SPREAD (INST -1 "f" "{domain_start}" "{domain_end}" "{deriv_bound}" "{max_right_clipped}" "x")
%|-           ((SPREAD (SPLIT -1)
%|-             ((THEN (LEMMA "{deriv_lemma}")
%|-               (SPREAD (INST -1 "f" "{domain_start}" "{domain_end}" "{deriv_bound}" "x" "{min_left}")
%|-                ((SPREAD (SPLIT -1)
%|-                  ((THEN (ASSERT) (EXPAND "f") (ASSERT)) (ASSERT) (PROPAX)
%|-                   (ASSERT) (PROPAX)))
%|-                 (THEN (EXPAND "ci") (ASSERT))
%|-                 (THEN (EXPAND "ci") (PROPAX)))))
%|-              (ASSERT) (PROPAX) (ASSERT) (PROPAX)))
%|-            (THEN (EXPAND "ci") (PROPAX))
%|-            (THEN (ASSERT) (EXPAND "ci") (PROPAX)))))
%|-         (THEN (ASSERT)
%|-          (SPREAD (CASE "{case_label4}") ((GRIND) (GRIND))))))))))))
"""


def generate_two_case_unifying_lemma_proof(
    domain_split,
    lemma_1,
    lemma_2,
    domain_type_1,
    domain_type_2,
    trajectory_function_1,
    trajectory_function_2,
):
    """
    Input Examples:
    - domain_split: 0
    - lemma_1: le_lo_case_0
    - lemma_2: ge_ro_case_1
    - domain_type_1: left_open(0)
    - domain_type_2: right_open(0)
    - trajectory_function_1: LAMBDA (x: real): COND x > 0 -> g(0), ELSE -> g(x) ENDCOND)
    - trajectory_function_2: LAMBDA (x: real): COND x >= 0 -> g(x), ELSE -> g(0) ENDCOND
    """

    return f"""%|- full_domain_soundness_lemma : PROOF
%|- (THEN (SKEEP) (SKOLETIN*) (FLATTEN)
%|-  (SPREAD (CASE "x<={domain_split}")
%|-   ((THEN (LEMMA "{lemma_1}") (INST -1 "xo" "yo" "g") (ASSERT) (EXPAND "f0")
%|-     (ASSERT)
%|-     (SPREAD
%|-      (CASE
%|-          "restrict[real, ({domain_type_1}), real](g) = (restrict[real, ({domain_type_1}), real]
%|-                 ({trajectory_function_1}))")
%|-      ((SPREAD (SPLIT -2)
%|-        ((PROPAX) (PROPAX) (PROPAX) (PROPAX) (PROPAX) (ASSERT) (ASSERT)
%|-         (THEN (INST 1 "x") (ASSERT))))
%|-       (THEN (DECOMPOSE-EQUALITY 1) (HIDE-ALL-BUT 1) (TYPEPRED "x!1")
%|-        (GRIND)))))
%|-    (THEN (LEMMA "{lemma_2}") (INST -1 "xo" "yo" "g") (ASSERT) (EXPAND "f1")
%|-     (SPREAD
%|-      (CASE "(restrict[real, ({domain_type_2}), real]
%|-                     ({trajectory_function_2})) = (restrict[real, ({domain_type_2}), real](g))")
%|-      ((SPREAD (SPLIT -2)
%|-        ((ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT)
%|-         (THEN (INST 1 "x") (ASSERT))))
%|-       (THEN (HIDE-ALL-BUT 1) (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!1")
%|-        (GRIND))))))))
%|- QED full_domain_soundness_lemma"""


def generate_three_case_unifying_lemma_proof(
    domain_split_1,
    domain_split_2,
    lemma_1,
    lemma_2,
    lemma_3,
    domain_type_1,
    domain_type_2,
    domain_type_3,
    trajectory_function_1,
    trajectory_function_2,
    trajectory_function_3,
):
    """
    Input Examples:
    - domain_split_1: 0
    - domain_split_2: 4
    - lemma_1: le_lo_case_0
    - lemma_2: ge_ro_case_1
    - domain_type_1: left_open(0)
    - domain_type_2: right_open(0)
    - trajectory_function_1: LAMBDA (x: real): COND x > 0 -> g(0), ELSE -> g(x) ENDCOND)
    - trajectory_function_2: LAMBDA (x: real): COND x >= 0 -> g(x), ELSE -> g(0) ENDCOND
    - trajectory_function_3: LAMBDA (x: real): COND x >= 0 -> g(x), ELSE -> g(0) ENDCOND
    """

    return f"""%|- full_domain_soundness_lemma : PROOF
%|- (THEN (SKEEP) (SKOLETIN*) (FLATTEN)
%|-  (SPREAD (CASE "x <= {domain_split_1}")
%|-   ((THEN (LEMMA "{lemma_1}") (INST -1 "xo" "yo" "g") (ASSERT) (EXPAND "f0")
%|-     (ASSERT)
%|-     (SPREAD
%|-      (CASE
%|-          "restrict[real, ({domain_type_1}), real](g) = (restrict[real, ({domain_type_1}), real]
%|-                 ({trajectory_function_1}))")
%|-      ((SPREAD (SPLIT -2)
%|-        ((PROPAX) (PROPAX) (PROPAX) (PROPAX) (PROPAX) (ASSERT) (ASSERT)
%|-         (THEN (INST 1 "x") (ASSERT))))
%|-       (THEN (DECOMPOSE-EQUALITY 1) (HIDE-ALL-BUT 1) (TYPEPRED "x!1")
%|-        (GRIND)))))
%|-    (SPREAD (CASE "x <= {domain_split_2}")
%|-     ((THEN (LEMMA "{lemma_2}") (INST -1 "xo" "yo" "g") (ASSERT)
%|-       (EXPAND "f1")
%|-       (SPREAD
%|-        (CASE
%|-            "restrict[real, ({domain_type_2}), real](g) = (restrict[real, ({domain_type_2}), real]({trajectory_function_2}))")
%|-        ((SPREAD (SPLIT -2)
%|-          ((ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT)
%|-           (ASSERT) (ASSERT) (THEN (INST 1 "x") (ASSERT))))
%|-         (THEN (ASSERT) (DECOMPOSE-EQUALITY 1) (HIDE-ALL-BUT 1) (TYPEPRED "x!1")
%|-          (GRIND))
%|-         (THEN (ASSERT) (HIDE-ALL-BUT 1) (GRIND)))))
%|-      (THEN (LEMMA "{lemma_3}") (INST -1 "xo" "yo" "g") (ASSERT)
%|-       (EXPAND "f2")
%|-       (SPREAD
%|-        (CASE
%|-            "restrict[real, ({domain_type_3}), real](g) = (restrict[real, ({domain_type_3}), real]({trajectory_function_3}))")
%|-        ((SPREAD (SPLIT -2)
%|-          ((ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT)
%|-           (THEN (ASSERT) (INST 1 "x") (ASSERT))))
%|-         (THEN (HIDE-ALL-BUT 1) (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!1")
%|-          (GRIND))))))))))
%|- QED full_domain_soundness_lemma
"""


def generate_tcc_postamble(case_tccs):
    all_tccs = "\n\n".join(case_tccs)
    return f"""
%|- *TCC* : PROOF (THEN (ASSERT) (GRIND)) QED

%|- mvt_gen_ge_TCC1 : PROOF (THEN (SKEEP) (LEMMA "connected_deriv_domain[(D)]") (ASSERT)) QED

%|- mvt_gen_ge_lo_TCC1 : PROOF (THEN (SKEEP) (LEMMA "left_open_dd") (INST -1 "C") (ASSERT) (GRIND)) QED

%|- mvt_gen_ge_lo_TCC2 : PROOF (THEN (SKEEP) (LEMMA "left_open_noe") (INST -1 "C") (ASSERT) (GRIND)) QED

%|- mvt_gen_ge_ro_TCC1 : PROOF (THEN (SKEEP) (LEMMA "right_open_dd") (INST -1 "C") (ASSERT) (GRIND)) QED

%|- mvt_gen_ge_ro_TCC2 : PROOF (THEN (SKEEP) (LEMMA "right_open_noe") (INST -1 "C") (ASSERT) (GRIND)) QED

%|- mvt_gen_ge_ci_TCC1 : PROOF (THEN (SKEEP) (LEMMA "ci_dd") (INST -1 "d1" "d2") (ASSERT) (GRIND)) QED

%|- mvt_gen_ge_ci_TCC2 : PROOF (THEN (SKEEP) (LEMMA "ci_noe") (INST -1 "d1" "d2") (ASSERT) (GRIND)) QED

{all_tccs}

end active_corner_certificate"""


# =============================================================================
# AUTOMATIC PROOF SCRIPT GENERATION
# =============================================================================
def generate_proof_calls(trajectory_expr, poly, domain, x=symbols("x"), y=symbols("y")):
    """
    Automatically generate proof script calls based on
    trajectory, polygon, and domain analysis.

    Args:
        trajectory_expr: Single-variable Sympy *expression* for trajectory. May be Piecewise.
        poly: Sympy Polygon object
        domain: Sympy Interval domain
        x, y: Sympy symbols for variables

    Returns:
        list: List of dictionaries containing parameters for unbounded_one_side_proof_script calls
    """
    # Trajectory is assumed to be function of x only
    func_var = x

    labels = list(string.ascii_lowercase[: len(poly.vertices)])
    labels = list(labels[1:] + labels[:1])
    verts, lines = verts_and_lines(labels, poly)

    # Generate premise and explicit conditions
    # Treat trajectory as a function f(x) rather than substituting the actual expression

    f = Function("f")(x)
    traj = y - f
    premise = generate_premise(lines, f)
    corner_pairs = generate_corner_pairs(labels)
    explicit = generate_explicit_disjunction(corner_pairs, traj, verts)

    # Compute polygon width using same method as construct_lemma
    max_offset = max([v.x for v in poly.vertices])
    min_offset = min([v.x for v in poly.vertices])

    # Calculate half-width (distance from centroid)
    half_width = (max_offset - min_offset) / 2

    # Get polygon angles and find transitions
    angles, _ = compute_polygon_angles(poly)
    set_of_transitions = set()

    if isinstance(trajectory_expr, Piecewise):
        # Handle piecewise trajectories
        for subtraj, subcond in trajectory_expr.as_expr_set_pairs():
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
            -y + trajectory_expr, angles, x, y, domain=domain
        )
        # Add left and right boundaries to check for notch there too
        left_bound = Point(domain.inf, trajectory_expr.subs(func_var, domain.inf))
        right_bound = Point(domain.sup, trajectory_expr.subs(func_var, domain.sup))

        set_of_transitions.update(subset_of_transitions)
        # only add if finite
        if left_bound.x.is_finite and left_bound.y.is_finite:
            set_of_transitions.add(left_bound)
        if right_bound.x.is_finite and right_bound.y.is_finite:
            set_of_transitions.add(right_bound)

    # Sort transitions by x coordinate and compute midpoints
    sorted_transitions = sorted(set_of_transitions, key=lambda point: point.x)
    func_var_transitions = [p.x for p in sorted_transitions]
    print(f"func_var_transitions: {func_var_transitions}")
    midpoints = np.convolve(func_var_transitions, [1, 1], "valid") / 2

    # Compute derivative and angles at midpoints
    if isinstance(trajectory_expr, Piecewise):
        # For piecewise, we need to handle each piece separately
        deriv_traj = trajectory_expr.diff(x)
    else:
        deriv_traj = diff(trajectory_expr, x)

    dydx_midpoints = []
    for val in midpoints:
        try:
            deriv_val = deriv_traj.subs(func_var, val)
            dydx_midpoints.append(deriv_val)
        except:
            # Handle cases where derivative evaluation fails
            print(f"Derivative evaluation failed at {val}")
            dydx_midpoints.append(0)  # Default to 0 if evaluation fails
    # subtract pi/2 to get angles in the range [-pi/2, pi/2]
    midpoint_angles = [atan2(d, 1) - pi / 2 for d in dydx_midpoints]

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
    # TODO maybe we add duplicates here? but maybe it doesn't matter.
    all_transition_points = func_var_transitions.copy()
    if domain_min.is_finite:
        all_transition_points.insert(0, domain_min)
    elif domain_min == -oo:
        all_transition_points.insert(0, -oo)
    if domain_max.is_finite:
        all_transition_points.append(domain_max)
    elif domain_max == oo:
        all_transition_points.append(oo)
    print(f"all_transition_points: {all_transition_points}")

    # Create intervals between all transition points
    var_intervals = list(zip(all_transition_points[:-1], all_transition_points[1:]))
    print(f"var_intervals: {var_intervals}")

    proof_calls = []
    deriv_statements = []
    notches = []
    trajectories = []
    num_cases = 0

    for i, var_interval in enumerate(var_intervals):
        interval_start, interval_end = var_interval

        # Skip if this interval is outside the domain
        if interval_end <= domain_min or interval_start >= domain_max:
            continue

        num_cases += 1

        # Determine if this interval is unbounded on left or right
        left_unbounded = not interval_start.is_finite
        right_unbounded = not interval_end.is_finite

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

        case_labels = []
        # Generate case label for unbounded domains
        if left_unbounded:  # Interval unbounded on left
            case_labels.append(f"xo + {half_width} <= {interval_end}")
        elif right_unbounded:  # right_unbounded: Interval unbounded on right
            case_labels.append(f"xo - {half_width} >= {interval_start}")
        else:
            # Closed interval case
            case_labels.append(
                f"xo - {half_width} >= {interval_start} AND xo + {half_width} <= {interval_end}"
            )
            case_labels.append(
                f"xo - {half_width} < {interval_start} AND xo + {half_width} <= {interval_end}"
            )
            case_labels.append(
                f"xo - {half_width} >= {interval_start} AND xo + {half_width} > {interval_end}"
            )
            case_labels.append(
                f"xo - {half_width} < {interval_start} AND xo + {half_width} > {interval_end}"
            )

        # Determine derivative lemma based on derivative sign at midpoint
        # For left_unbounded domains: use left_open lemmas
        # For right_unbounded domains: use right_open lemmas
        deriv_sign = "ge" if deriv_at_midpoint >= 0 else "le"
        if left_unbounded:
            domain_type = "lo"
        elif right_unbounded:
            domain_type = "ro"
        else:
            domain_type = "ci"
        deriv_lemma = f"mvt_gen_{deriv_sign}_{domain_type}"

        # Polygon bounds
        max_right = f"xo + {half_width}"
        min_left = f"xo - {half_width}"

        # bounds on "xo" for the case where we clip the trajectory
        if left_unbounded:
            min_left_clipped = min_left
            max_right_clipped = interval_end
        elif right_unbounded:
            min_left_clipped = interval_start
            max_right_clipped = max_right
        else:  # closed interval
            min_left_clipped = interval_start
            max_right_clipped = interval_end

        # Compute domain bounds for lemma generation
        if left_unbounded:
            domain_start = None
            domain_end = interval_end
            notch = premise.subs(x, interval_end)
            domain_definition = f"left_open({domain_end})"
        elif right_unbounded:
            domain_start = interval_start
            domain_end = None
            notch = premise.subs(x, interval_start)
            domain_definition = f"right_open({domain_start})"
        else:
            domain_start = interval_start
            domain_end = interval_end
            notch = premise.subs(x, interval_start) | premise.subs(x, interval_end)
            domain_definition = f"ci({domain_start}, {domain_end})"

        # Determine deriv_clause1 based on derivative sign
        if left_unbounded or right_unbounded:
            endpoint = interval_end if left_unbounded else interval_start
            deriv_at_endpoint = deriv_traj.subs(x, endpoint)
            deriv_direction = ">=" if deriv_at_midpoint >= deriv_at_endpoint else "<="
            deriv_clause1 = f"{deriv_direction} {deriv_at_endpoint}"
            deriv_clause2 = None
            deriv_bound1 = deriv_at_endpoint
            deriv_bound2 = None
        else:
            # TODO see if this still works for intervals with piecewise behavior i.e.
            # parabola-then-straight
            deriv_at_beginning = deriv_traj.subs(x, interval_start)
            deriv_at_end = deriv_traj.subs(x, interval_end)
            deriv_direction_1 = (
                ">=" if deriv_at_beginning <= deriv_at_midpoint else "<="
            )
            deriv_direction_2 = "<=" if deriv_direction_1 == ">=" else ">="
            deriv_clause1 = f"{deriv_direction_1} {deriv_at_beginning}"
            deriv_clause2 = f"{deriv_direction_2} {deriv_at_end}"
            deriv_bound1 = deriv_at_beginning
            deriv_bound2 = deriv_at_end
        print(f"deriv_clause1: {deriv_clause1}")
        print(f"deriv_clause2: {deriv_clause2}")

        # Generate lemma name for this interval
        lemma_name = f"{deriv_sign}_{domain_type}_case_{i}"

        # Generate the lemma
        lemma_text, deriv_list, clipped_trajectory = construct_lemma(
            verts,
            premise,
            explicit,
            lemma_name,
            deriv_clause1=deriv_clause1,
            deriv_clause2=deriv_clause2,
            domain_start=domain_start,
            domain_end=domain_end,
            notch=notch,
        )

        proof_call = {
            "lemma_name": lemma_name,
            "case_labels": case_labels,
            "deriv_lemma": deriv_lemma,
            "max_right": max_right,
            "min_left": min_left,
            "domain_definition": domain_definition,
            "domain_start": domain_start,
            "domain_end": domain_end,
            "max_right_clipped": max_right_clipped,
            "min_left_clipped": min_left_clipped,
            "interval": var_interval,
            "active_corner": active_corner_offset,
            "lemma_text": lemma_text,
            "deriv_bound1": deriv_bound1,
            "deriv_bound2": deriv_bound2,
        }

        proof_calls.append(proof_call)
        deriv_statements.extend(deriv_list)
        notches.append(notch)
        trajectories.append(clean_trajectory(clipped_trajectory))

    domain_splits = []
    for i in range(len(var_intervals) - 1):
        if var_intervals[i][1] == var_intervals[i + 1][0]:
            domain_splits.append(var_intervals[i][1])

    # Only generate unifying lemma if we can handle creating the proof too
    unifying_lemma_statement = generate_unifying_lemma(
        verts, premise, explicit, deriv_statements, notches, trajectories
    )
    if num_cases == 2:
        lemma_1 = proof_calls[0]["lemma_name"]
        lemma_2 = proof_calls[1]["lemma_name"]
        domain_split = domain_splits[0]
        domain_type_1 = proof_calls[0]["domain_definition"]
        domain_type_2 = proof_calls[1]["domain_definition"]
        trajectory_function_1 = trajectories[0]
        trajectory_function_2 = trajectories[1]
        unifying_proof = generate_two_case_unifying_lemma_proof(
            domain_split,
            lemma_1,
            lemma_2,
            domain_type_1,
            domain_type_2,
            trajectory_function_1,
            trajectory_function_2,
        )
        unifying_lemma_and_proof = unifying_lemma_statement + "\n\n" + unifying_proof
    elif num_cases == 3:
        lemma_1 = proof_calls[0]["lemma_name"]
        lemma_2 = proof_calls[1]["lemma_name"]
        lemma_3 = proof_calls[2]["lemma_name"]
        domain_split_1 = domain_splits[0]
        domain_split_2 = domain_splits[1]
        domain_type_1 = proof_calls[0]["domain_definition"]
        domain_type_2 = proof_calls[1]["domain_definition"]
        domain_type_3 = proof_calls[2]["domain_definition"]
        trajectory_function_1 = trajectories[0]
        trajectory_function_2 = trajectories[1]
        trajectory_function_3 = trajectories[2]
        unifying_proof = generate_three_case_unifying_lemma_proof(
            domain_split_1,
            domain_split_2,
            lemma_1,
            lemma_2,
            lemma_3,
            domain_type_1,
            domain_type_2,
            domain_type_3,
            trajectory_function_1,
            trajectory_function_2,
            trajectory_function_3,
        )
        unifying_lemma_and_proof = unifying_lemma_statement + "\n\n" + unifying_proof
    else:
        unifying_lemma_and_proof = unifying_lemma_statement

    return proof_calls, unifying_lemma_and_proof


def generate_lemmas_from_calls(proof_calls):
    """
    Generate lemma strings from proof call parameters.

    Args:
        proof_calls: List of dictionaries from generate_unbounded_proof_calls

    Returns:
        list: List of lemma strings
    """
    lemmas = []

    for call_params in proof_calls:
        lemma_text = call_params.get("lemma_text", "")
        if lemma_text:
            lemmas.append(lemma_text)

    return lemmas


def generate_tcc_from_call(proof_call):
    """
    Generate a TCC from a proof call.
    """
    tcc_name = f"{proof_call['lemma_name']}_TCC"

    if proof_call["domain_definition"].startswith("left_open"):
        assert proof_call["domain_start"] is None
        lemma_name = "left_open"
        inst_call = f"(INST -1 \"{proof_call['domain_end']}\")"
    elif proof_call["domain_definition"].startswith("right_open"):
        assert proof_call["domain_end"] is None
        lemma_name = "right_open"
        inst_call = f"(INST -1 \"{proof_call['domain_start']}\")"
    elif proof_call["domain_definition"].startswith("ci"):
        lemma_name = "ci"
        inst_call = (
            f"(INST -1 \"{proof_call['domain_start']}\" \"{proof_call['domain_end']}\")"
        )
    else:
        raise ValueError(
            f"Unknown domain definition: {proof_call['domain_definition']}"
        )

    return f"""%|- {tcc_name}3 : PROOF (THEN (SKEEP) (LEMMA "{lemma_name}_dd") {inst_call} (ASSERT) (GRIND)) QED

%|- {tcc_name}4 : PROOF (THEN (SKEEP) (LEMMA "{lemma_name}_noe") {inst_call} (ASSERT) (GRIND)) QED"""


def generate_corner_pairs(labels):
    """
    Generate corner pairs for a regular polygon.

    Only supports even-numbered polygons since only they have true opposite vertices.
    For odd-numbered polygons, returns an empty list.

    Args:
        labels: List of vertex labels (e.g., ['a', 'b', 'c', 'd'])

    Returns:
        list: List of tuples representing opposite vertex pairs, or empty list for odd polygons
    """
    n = len(labels)
    # Only generate corner pairs for even-numbered polygons
    # Odd-numbered polygons don't have true opposite vertices
    if n < 4 or n % 2 != 0:
        return []

    pairs = []
    # For even n, vertices are exactly opposite
    offset = n // 2

    for i in range(n):
        opposite_idx = (i + offset) % n
        # Only add each pair once (avoid duplicates)
        if i < opposite_idx:
            pairs.append((labels[i], labels[opposite_idx]))

    return pairs


def generate_complete_proof_package(trajectory, poly, domain, lemma_name="testlemma"):
    """
    Generate a complete proof package including lemmas and proof scripts.

    This function demonstrates how to use the new lemma generation functionality
    similar to the notebook example pattern.

    Args:
        trajectory: The trajectory expression (e.g., Function('f')(x))
        poly: The polygon object
        domain: The domain intervals
        lemma_name: Base name for the lemma

    Returns:
        dict: Dictionary containing lemmas and proof scripts
    """

    # Create trajectory expression similar to notebook example
    x, y = symbols("x y")
    if isinstance(trajectory, str):
        # If trajectory is a string, assume it's a function name
        traj_expr = Function(trajectory)(x)
    else:
        traj_expr = trajectory

    # Generate proof calls
    proof_calls, unifying_lemma = generate_proof_calls(traj_expr, poly, domain, x, y)

    # Generate lemmas
    lemmas = generate_lemmas_from_calls(proof_calls)

    # Generate proof scripts
    proof_scripts = generate_proof_scripts_from_calls(proof_calls)

    # Generate TCCs
    tccs = [generate_tcc_from_call(call) for call in proof_calls]

    # Create a complete package
    package = {
        "proof_calls": proof_calls,
        "lemmas": lemmas,
        "proof_scripts": proof_scripts,
        "trajectory": traj_expr,
        "polygon": poly,
        "domain": domain,
        "tccs": tccs,
        "unifying_lemma": unifying_lemma,
    }

    return package


def print_prooflite(package) -> str:
    """
    Print the proof package in a readable format.
    """
    lemmas = package["lemmas"]
    proof_scripts = package["proof_scripts"]
    lemma_names = [lemma.split(":")[0] for lemma in lemmas]

    with open("certificate_preamble.txt") as file:
        s = file.read()
    for i in range(len(lemmas)):
        s += lemmas[i] + "\n"
        s += f"%|- {lemma_names[i]} : PROOF\n"
        s += proof_scripts[i]
        s += f"%|- QED {lemma_names[i]}\n\n\n"

    # only print unifying lemma if we can handle creating the proof too
    if package["unifying_lemma"] is not None:
        s += package["unifying_lemma"] + "\n\n"

    s += generate_tcc_postamble(package["tccs"])

    return s.strip()


def log_proof_to_file(package, filename):
    """
    Log the proof package to a file.
    """
    with open(filename, "w") as f:
        f.write(print_prooflite(package))


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
        if call_params["domain_definition"].startswith("ci"):
            assert len(call_params["case_labels"]) == 4
            script = bounded_proof_script(
                case_label1=call_params["case_labels"][0],
                case_label2=call_params["case_labels"][1],
                case_label3=call_params["case_labels"][2],
                case_label4=call_params["case_labels"][3],
                deriv_lemma=call_params["deriv_lemma"],
                max_right=call_params["max_right"],
                min_left=call_params["min_left"],
                max_right_clipped=call_params["max_right_clipped"],
                min_left_clipped=call_params["min_left_clipped"],
                domain_start=call_params["domain_start"],
                domain_end=call_params["domain_end"],
                deriv_bound=call_params["deriv_bound1"],
            )
        else:
            assert len(call_params["case_labels"]) == 1

            if call_params["domain_definition"].startswith("left_open"):
                domain_bound = call_params["domain_end"]
            elif call_params["domain_definition"].startswith("right_open"):
                domain_bound = call_params["domain_start"]
            else:
                raise ValueError(
                    f"Unknown domain definition: {call_params['domain_definition']}"
                )
            assert domain_bound is not None, f"domain_bound is None for {call_params}"

            script = unbounded_one_side_proof_script(
                case_label=call_params["case_labels"][0],
                deriv_lemma=call_params["deriv_lemma"],
                max_right=call_params["max_right"],
                min_left=call_params["min_left"],
                domain_definition=call_params["domain_definition"],
                max_right_clipped=call_params["max_right_clipped"],
                min_left_clipped=call_params["min_left_clipped"],
                domain_bound=domain_bound,
                deriv_bound=call_params["deriv_bound1"],
            )
        proof_scripts.append(script)

    return proof_scripts
