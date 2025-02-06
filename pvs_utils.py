import numpy as np
import matplotlib.pyplot as plt
from sympy import expand, Line, Polygon, symbols, Point, Function


def plot_polygon(poly: Polygon, labels=[]):
    # Draw a polygon by plotting vertices as points and edges as lines
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


def generate_premise(lines: dict, traj_expr):
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
    explicit = False
    for c1, c2 in corner_pairs:
        corner_traj_1 = expand(traj.subs(y, yo - verts[c1].y).subs(x, xo - verts[c1].x))
        corner_traj_2 = expand(traj.subs(y, yo - verts[c2].y).subs(x, xo - verts[c2].x))
        pair_clause = ((corner_traj_1 <= 0) & (corner_traj_2 >= 0)) | (
            (corner_traj_1 >= 0) & (corner_traj_2 <= 0)
        )
        explicit |= pair_clause
    return explicit


def sympy_to_pvs(clause: str):
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


def build_line(x, y, coefs):
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
    if type(expr) == dict:
        return make_relative_dict(expr, traj_expr)
    if traj_expr is None:
        return expr.subs(x, xo - x).subs(y, yo - alpha * x)
    else:
        return expr.subs(x, xo - x).subs(y, yo - traj_expr)


def make_relative_dict(d, te):
    return {k: make_relative(v, te) for k, v in d.items()}


def verts_and_lines(
    vert_names: list[str], poly: Polygon, x=symbols("x"), y=symbols("y")
):
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


# TODO(nishant): automatically extract exists-clause bounds from the premise
def construct_lemma(
    active_exists_premise,
    active_corner_condition,
    lemma_name,
    deriv_clause1="derivable?[real](f)",
    deriv_clause2=" >= 0)",
    deriv_clause3=None,
    exists_upper=None,
    exists_lower=None,
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
    exists_clause = f"(EXISTS (x : real) :\n    ({sympy_to_pvs(str(active_exists_premise))}) AND\n    {exists_upper} >= x AND {exists_lower} <= x)"

    if deriv_clause3:
        full_preamble = f"{deriv_clause1} AND\n    (FORALL(x:real): deriv[real](f)(x) {deriv_clause2}) AND\n    (FORALL(x:real): deriv[real](f)(x) {deriv_clause3}) AND\n    {exists_clause}"
    else:
        full_preamble = f"{deriv_clause1} AND\n    (FORALL(x:real): deriv[real](f)(x) {deriv_clause2}) AND\n    {exists_clause}"

    lemma_str = f"""{lemma_name}: LEMMA
    FORALL(f:[real-> real],xo,yo:real):
    {full_preamble}
    IMPLIES
    {sympy_to_pvs(str(active_corner_condition))}
"""
    return lemma_str


class ProofNode:
    def __init__(self, command, args=None):
        self.command = command
        self.args = args if args else []
        self.children = []

    def __repr__(self):
        args_str = f" {self.args}" if self.args else ""
        children_str = f" children={len(self.children)}" if self.children else ""
        return f"ProofNode({self.command}{args_str}{children_str})"

    def add_child(self, node):
        self.children.append(node)
        return node

    def add_children(self, nodes):
        self.children.extend(nodes)
        return nodes

    def generate(self, indent="%|- ", depth=0, n_spaces=4):
        if self.command == "()":  # Special handling for wrapper nodes
            # Generate children with extra parentheses
            child_proofs = [c.generate(indent, depth) for c in self.children]
            return f"({' '.join(child_proofs)})"
        base_indent = indent + " " * depth * n_spaces

        # Base case - no children
        if not self.children:
            return f"({self.command}{' ' + ' '.join(str(arg) for arg in self.args) if self.args else ''})"

        # Multiple children or nested structures
        result = f"({self.command}{' ' + ' '.join(str(arg) for arg in self.args) if self.args else ''}"

        child_proofs = [c.generate(indent, depth + 1) for c in self.children]
        child_str = f"\n{base_indent} ".join(child_proofs)

        return f"{result}\n{base_indent} {child_str})"


class ProofBuilder:
    def __init__(self):
        self.indent = "%|- "

    def __repr__(self):
        return f"ProofBuilder(indent='{self.indent}')"

    def create_then_sequence(self, *steps):
        then = ProofNode("THEN")
        for step in steps:
            if isinstance(step, list):
                for s in step:
                    then.add_child(s)
            else:
                then.add_child(step)
        return then

    def create_spread_split(self, branches):
        spread = ProofNode("SPREAD")
        split = ProofNode("SPLIT", ["-1"])
        spread.add_child(split)

        # Create a wrapper node to enclose all branch nodes
        wrapper = ProofNode("()")
        for branch in branches:
            wrapper.add_child(branch[0])
        spread.add_child(wrapper)

        return spread

    def create_mvt_step(self, mvt_lemma):
        return [
            ProofNode("LEMMA", [f'"{mvt_lemma}"']),
            ProofNode("INST?"),
            ProofNode("ASSERT"),
        ]

    def create_inst_step(self, term1, term2):
        return [
            ProofNode("INST", ["-1", f'"{term1}"', f'"{term2}"']),
            ProofNode("ASSERT"),
        ]

    def build_rect_proof(self, lemma_name, mvt_lemma, inst_terms):
        script = ProofScript(lemma_name)

        root = self.create_then_sequence(
            ProofNode("SKEEP*"), *self.create_mvt_step(mvt_lemma)
        )

        # Inner branch with THEN sequence
        inner_then = self.create_then_sequence(
            *self.create_inst_step(inst_terms[1], inst_terms[2])
        )

        # Inner SPREAD with branches in list
        inner_spread = self.create_spread_split([[inner_then], [ProofNode("PROPAX")]])

        # Main branch with THEN sequence
        main_branch = self.create_then_sequence(
            ProofNode("ASSERT"),
            *self.create_inst_step(inst_terms[0], inst_terms[1]),
            *self.create_mvt_step(mvt_lemma),
            inner_spread,
        )

        # Outer SPREAD with branches in list
        outer_spread = self.create_spread_split([[main_branch], [ProofNode("PROPAX")]])
        root.add_child(outer_spread)
        script.root = root
        return script


class ProofScript:
    def __init__(self, lemma_name):
        self.lemma_name = lemma_name
        self.root = None

    def generate(self):
        try:
            lines = []
            lines.append(f"%|- {self.lemma_name} : PROOF")
            lines.append(f"%|- {self.root.generate()}")
            lines.append(f"%|- QED {self.lemma_name}")
            return "\n".join(lines)
        except:
            print(self.lemma_name)
            print(self.root)
