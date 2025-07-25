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


class ProofNode:
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
        self.children.append(node)
        return node

    def add_children(self, nodes):
        self.children.extend(nodes)
        return nodes

    def generate(self, indent="%|- ", depth=0, n_spaces=4):
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

    def create_spread_case(self, case_str, branches):
        spread = ProofNode("SPREAD")
        case = ProofNode("CASE", [f'"{case_str}"'])
        spread.add_child(case)
        for branch in branches:
            spread.add_child(branch)
        return spread

    def create_spread_split(self, branches):
        spread = ProofNode("SPREAD")
        split = ProofNode("SPLIT", ["-1"])
        spread.add_child(split)
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

    def create_inst_step(self, *args, **kwargs):
        # args: positional arguments for INST, kwargs: keyword arguments (e.g., :WHERE 1)
        return [ProofNode("INST", list(args), kwargs), ProofNode("ASSERT")]

    def create_instq_step(self, *args, **kwargs):
        # For INST? with optional :WHERE, etc.
        return [ProofNode("INST?", list(args), kwargs), ProofNode("ASSERT")]

    def build_example_proof(self):
        # Example: builds a proof tree matching the structure in bound22_rect_function_bounded
        # (THEN (SKEEP*) (SKOLETIN*) (FLATTEN) (SKEEP) (ASSERT)
        #   (SPREAD (CASE "xo-2 >=0")
        #     (...)
        #     (...)
        #   )
        # )
        then = self.create_then_sequence(
            ProofNode("SKEEP*"),
            ProofNode("SKOLETIN*"),
            ProofNode("FLATTEN"),
            ProofNode("SKEEP"),
            ProofNode("ASSERT"),
        )
        # Branch 1 (abbreviated)
        branch1 = self.create_then_sequence(
            ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
            ProofNode("INST", ["-1", '"f"', '"xo + 2"', '"x"', '"0"']),
            ProofNode(
                "SPREAD",
                children=[
                    ProofNode("SPLIT", ["-1"]),
                    ProofNode(
                        "()",
                        children=[
                            self.create_then_sequence(
                                ProofNode("ASSERT"),
                                ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
                                ProofNode(
                                    "INST", ["-1", '"f"', '"x"', '"xo-2"', '"0"']
                                ),
                                ProofNode("ASSERT"),
                                ProofNode("SKEEP"),
                                ProofNode("INST?"),
                                ProofNode("ASSERT"),
                            ),
                            ProofNode("PROPAX"),
                            ProofNode("ASSERT"),
                            self.create_then_sequence(
                                ProofNode("SKEEP"),
                                ProofNode("INST?", [], {"WHERE": 1}),
                                ProofNode("ASSERT"),
                            ),
                        ],
                    ),
                ],
            ),
        )
        # Branch 2 (abbreviated)
        branch2 = self.create_then_sequence(
            ProofNode("ASSERT"),
            ProofNode("EXPAND", ['"f"']),
            ProofNode("ASSERT"),
            ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
            ProofNode("INST", ["-1", '"f"', '"xo+2"', '"x"', '"0"']),
            ProofNode(
                "SPREAD",
                children=[
                    ProofNode("SPLIT", ["-1"]),
                    ProofNode(
                        "()",
                        children=[
                            self.create_then_sequence(
                                ProofNode("ASSERT"),
                                ProofNode("EXPAND", ['"f"']),
                                ProofNode(
                                    "SPREAD",
                                    children=[
                                        ProofNode("CASE", ['"x=0"']),
                                        self.create_then_sequence(
                                            ProofNode("ASSERT"),
                                            ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
                                            ProofNode(
                                                "INST",
                                                ["-1", '"f"', '"x"', '"0"', '"0"'],
                                            ),
                                            ProofNode("ASSERT"),
                                            ProofNode("SKEEP"),
                                            ProofNode("INST?"),
                                            ProofNode("ASSERT"),
                                        ),
                                    ],
                                ),
                            ),
                            ProofNode("EXPAND", ['"f"', "1"]),
                            ProofNode("PROPAX"),
                            ProofNode("ASSERT"),
                            self.create_then_sequence(
                                ProofNode("SKEEP"),
                                ProofNode("INST?", [], {"WHERE": 1}),
                                ProofNode("ASSERT"),
                            ),
                        ],
                    ),
                ],
            ),
        )
        spread = self.create_spread_case("xo-2 >=0", [branch1, branch2])
        then.add_child(spread)
        return then

    def build_rect_proof(
        self,
        lemma_name,
        mvt_lemma,
        inst_main,
        spread_cases,
        preamble_nodes=None,
    ):
        """
        Flexible rectangle proof builder (structure-driven).
        lemma_name: name of the lemma
        mvt_lemma: name of the MVT lemma to use (e.g., "mvt_gen_ge_bound")
        inst_main: list of terms for the main INST (e.g., ["f", "xo + 2", "x", "0"])
        spread_cases: list of dicts, each with:
            - 'case_label': label for the CASE
            - 'split_branches': list of dicts, each with 'type' and optional params
        preamble_nodes: optional list of ProofNode for THEN preamble (default: standard)
        Returns: ProofScript
        """

        def make_branch(branch):
            btype = branch.get("type")
            if btype == "then":
                steps = branch.get("steps", [])
                nodes = []
                for step in steps:
                    cmd = step["cmd"]
                    args = step.get("args", [])
                    kwargs = step.get("kwargs", {})
                    nodes.append(ProofNode(cmd, args, kwargs))
                return self.create_then_sequence(*nodes)
            elif btype == "assert-then":
                steps = branch.get("steps", [])
                nodes = [ProofNode("ASSERT")]
                for step in steps:
                    cmd = step["cmd"]
                    args = step.get("args", [])
                    kwargs = step.get("kwargs", {})
                    nodes.append(ProofNode(cmd, args, kwargs))
                return self.create_then_sequence(*nodes)
            elif btype == "propax":
                return ProofNode("PROPAX")
            elif btype == "assert":
                return ProofNode("ASSERT")
            else:
                raise ValueError(f"Unknown branch type: {btype}")

        script = ProofScript(lemma_name)
        if preamble_nodes is None:
            preamble_nodes = [
                ProofNode("SKEEP*"),
                ProofNode("SKOLETIN*"),
                ProofNode("FLATTEN"),
                ProofNode("SKEEP"),
            ]
        root = self.create_then_sequence(
            *preamble_nodes,
            ProofNode(
                "SPREAD",
                children=[
                    ProofNode("CASE", [f'"{case['case_label']}"']) if i == 0 else None
                    for i, case in enumerate(spread_cases)
                ]
                + [
                    ProofNode(
                        "()",
                        children=[
                            self.create_then_sequence(
                                ProofNode("LEMMA", [f'"{mvt_lemma}"']),
                                ProofNode(
                                    "INST", ["-1"] + [f'"{t}"' for t in inst_main]
                                ),
                            ),
                            ProofNode(
                                "SPREAD",
                                children=[
                                    ProofNode("SPLIT", ["-1"]),
                                    ProofNode(
                                        "()",
                                        children=[
                                            *(
                                                make_branch(branch)
                                                for branch in case["split_branches"]
                                            )
                                        ],
                                    ),
                                ],
                            ),
                        ],
                    )
                    for case in spread_cases
                ],
            ),
        )
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
        except Exception as e:
            print(self.lemma_name)
            print(self.root)
            print(e)
