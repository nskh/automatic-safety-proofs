"""
Modular implementation of unifying proof generation.

This module provides a flexible interface for generating PVS proof scripts
that can handle an arbitrary number of cases, replacing the fixed-case
implementations in pvs_utils.py.
"""

from typing import List, Dict, Optional, Union


def format_pvs_proof_line(line: str, max_line_length: int = 100) -> str:
    """
    Format a long PVS proof line by breaking it into multiple lines at logical points.

    Breaks lines after closing parentheses when followed by opening parentheses,
    or when lines exceed max_line_length.

    Args:
        line: The PVS proof line to format (without %|- prefix)
        max_line_length: Maximum line length before breaking

    Returns:
        Formatted multi-line string with proper indentation
    """
    # If line is short enough, return as-is
    if len(line) <= max_line_length:
        return line

    # Remove leading/trailing whitespace
    line = line.strip()

    # If it's a simple command, return as-is
    if not line.startswith("("):
        return line

    # Break at logical points: after "))" sequences, before nested (SPREAD, (CASE, etc.
    # Track string literals to avoid breaking inside them
    result_lines = []
    current_line = ""
    i = 0
    paren_depth = 0
    in_string = False
    escape_next = False

    while i < len(line):
        char = line[i]
        current_line += char

        # Track string literals (don't break inside them)
        if escape_next:
            escape_next = False
        elif char == '"':
            in_string = not in_string
        elif char == "\\" and in_string:
            escape_next = True

        # Only consider breaking if we're not inside a string
        if not in_string:
            if char == "(":
                paren_depth += 1
                # Check if we should break before this opening paren when line is getting long
                # Look back to see if we have a pattern like ") (" or space before "("
                if len(current_line) > max_line_length * 0.6:
                    # Check if previous char is ")" or space, indicating a good break point
                    if len(current_line) > 1:
                        prev_char = current_line[-2] if len(current_line) >= 2 else ""
                        if prev_char in (")", " "):
                            # Move the "(" to next line
                            result_lines.append(current_line[:-1].rstrip())
                            current_line = "   " + char
            elif char == ")":
                paren_depth -= 1

                # Look ahead to see what comes next
                j = i + 1
                # Skip whitespace
                while j < len(line) and line[j] == " ":
                    j += 1

                if j < len(line):
                    next_char = line[j]
                    # Break after "))" if followed by "(" or if line is getting long
                    should_break = False
                    if i > 0 and line[i - 1] == ")":
                        # We have "))" - break if followed by "("
                        if next_char == "(":
                            should_break = True
                    elif len(current_line) > max_line_length * 0.7:
                        # Line is getting long - break if next is opening paren
                        if next_char == "(":
                            should_break = True

                    if should_break:
                        result_lines.append(current_line.rstrip())
                        current_line = "   "  # Indent continuation
                        i = j - 1  # Will increment at end of loop
                        continue

        i += 1

    # Add remaining line
    if current_line.strip():
        result_lines.append(current_line.rstrip())

    # Join with newlines and proper prefix
    if len(result_lines) == 1:
        return result_lines[0]
    else:
        return "\n%|-    ".join(result_lines)


def generate_one_case_unifying_proof(lemma_name: str):
    """
    Generate proof for a single case.

    Args:
        lemma_name: Name of the lemma to use

    Returns:
        PVS proof script as string
    """
    return f"""%|- full_domain_soundness_lemma : PROOF
%|- (THEN (SKEEP) (SKOLETIN*) (FLATTEN) (EXPAND "g_1") (LEMMA "{lemma_name}")
%|-  (ASSERT))
%|- QED full_domain_soundness_lemma
"""


def generate_unifying_proof(
    cases: List[Dict[str, str]],
    full_traj: str,
    piecewise_splits_map: Optional[Dict[int, List[str]]] = None,
    use_case_statements: bool = True,
) -> str:
    """
    Generate a unifying proof for an arbitrary number of cases.

    This is a modular replacement for generate_two_case_unifying_proof and
    generate_three_case_unifying_proof that supports any number of cases.

    Args:
        cases: List of case dictionaries, each containing:
            - domain: Domain definition (e.g., "left_open(0)")
            - trajectory: Trajectory piece for this domain (e.g., "LAMBDA (x: real): ...")
            - domain_type: Domain type name (e.g., "left_open")
            - domain_split: Domain split point (e.g., "0") - optional, used for domain_dd lemma
        full_traj: Full trajectory expression
        piecewise_splits_map: Optional dict mapping case index to list of piecewise splits
            for that case.

            **Piecewise Splits Interface:**
            The piecewise_splits_map allows you to specify which domains have piecewise
            boundary points (discontinuities or changes in the trajectory function).

            Format: {domain_index: [split_value1, split_value2, ...]}

            Example:
            - If you have 3 domains (indices 0, 1, 2)
            - Domain 2 (the last one) has piecewise boundaries at x=4 and x=8
            - You would pass: piecewise_splits_map = {2: ["4", "8"]}

            Note: For 3-case proofs, piecewise splits are only used for the last domain
            (matching original behavior). For 4+ cases, splits can be specified for any domain.

            If None, will attempt to infer from cases (for backward compatibility).
        use_case_statements: Whether to use CASE statements for trajectory equality
            (True for 3+ cases, False for 2 cases)

    Returns:
        PVS proof script as string

    Example:
        cases = [
            {
                "domain": "left_open(0)",
                "trajectory": "LAMBDA (x: real): COND x > 0 -> g(0), ELSE -> g(x) ENDCOND",
                "domain_type": "left_open",
                "domain_split": "0"
            },
            {
                "domain": "right_open(0)",
                "trajectory": "LAMBDA (x: real): COND x >= 0 -> g(x), ELSE -> g(0) ENDCOND",
                "domain_type": "right_open",
                "domain_split": "0"
            }
        ]
        proof = generate_unifying_proof(cases, full_traj="...")
    """
    num_cases = len(cases)

    if num_cases == 0:
        return "% No cases provided"

    if num_cases == 1:
        # For single case, use simpler proof structure
        return "% Single case - use generate_one_case_unifying_proof instead"

    # Extract domain types and prepare trajectory equality statements
    domain_types = [case["domain_type"] for case in cases]
    domains = [case["domain"] for case in cases]
    trajectories = [case["trajectory"] for case in cases]

    # Generate PVS trajectory equality statements for CASE statements
    pvs_traj_equalities = []
    for i, (domain, traj) in enumerate(zip(domains, trajectories)):
        pvs_traj_eq = (
            f"(LAMBDA(s: ({domain})): {full_traj.replace('x', 's')}) = "
            f"(LAMBDA (s: ({domain})): {traj.replace('x', 's')})"
        )
        pvs_traj_equalities.append(pvs_traj_eq)

    # Determine number of PROPAX calls at the start
    # Pattern: 10 for 2 cases, 16 for 3 cases
    # For 4+ cases, we'll use 4*num_cases + 2 as a reasonable default
    if num_cases == 2:
        num_propax = 10
    elif num_cases == 3:
        num_propax = 16
    else:
        # For 4+ cases, estimate based on pattern
        num_propax = 4 * num_cases + 2

    # Generate proof branches for each case
    proof_branches = []

    if num_cases == 2:
        # Two-case proof structure (simpler, no CASE statements)
        domain_split = cases[0].get("domain_split", "")
        domain_type_1 = domain_types[0]
        domain_type_2 = domain_types[1]

        proof_branches = [
            f'(THEN (ASSERT) (HIDE-ALL-BUT (-7 1)) (REPLACE -1) (EXPAND "restrict")',
            f'     (WITH-TCCS (DERIVABLE 1)) (LEMMA "{domain_type_1}_dd") (INST -1 "{domain_split}"))',
            f'(THEN (HIDE-ALL-BUT (-7 1)) (REPLACE -1) (EXPAND "restrict") (SKEEP)',
            f"     (SPREAD (DERIV)",
            f'      ((THEN (TYPEPRED "x!1") (EXPAND "{domain_type_1}") (ASSERT))',
            f'       (THEN (LEMMA "{domain_type_1}_dd") (INST -1 "{domain_split}")))))',
            f'(THEN (HIDE-ALL-BUT (-7 1)) (REPLACE -1) (EXPAND "restrict") (DERIVABLE))',
            f'(THEN (HIDE-ALL-BUT (-7 1)) (REPLACE -1) (EXPAND "restrict") (SKEEP) (DERIV)',
            f'     (TYPEPRED "x!1") (EXPAND "{domain_type_2}") (ASSERT))',
        ]
    else:
        # Three or more cases: use CASE statements for trajectory equality
        # Check if we have piecewise splits (only used for last domain in 3-case)
        piecewise_splits = []
        if piecewise_splits_map and (num_cases - 1) in piecewise_splits_map:
            piecewise_splits = piecewise_splits_map[num_cases - 1]

        # For 3-case proofs, the original code returns early if no piecewise splits
        # This is a limitation we preserve for backward compatibility
        if num_cases == 3 and not piecewise_splits:
            return "% skipping unifying lemma proof."

        for i, (domain_type, pvs_traj_eq) in enumerate(
            zip(domain_types, pvs_traj_equalities)
        ):
            domain_split = cases[i].get("domain_split", "")

            # Generate branches for this domain
            if i == 0:
                # First domain: DERIVABLE branch with CASE
                # Match original formatting with explicit line breaks
                proof_branches.append(
                    f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict")\n'
                    f"     (SPREAD (DERIVABLE)\n"
                    f"      ((SPREAD\n"
                    f'        (CASE "{pvs_traj_eq}")\n'
                    f'        ((THEN (REPLACE -1) (HIDE -1) (DERIVABLE) (LEMMA "{domain_type}_dd")\n'
                    f'          (INST -1 "{domain_split}"))\n'
                    f'         (THEN (DECOMPOSE-EQUALITY 1) (HIDE 2) (TYPEPRED "x!1")\n'
                    f'          (EXPAND "{domain_type}" -1) (ASSERT))))\n'
                    f'       (THEN (LEMMA "{domain_type}_dd") (INST -1 "{domain_split}")))))'
                )

                # First domain: DERIV branch with CASE
                proof_branches.append(
                    f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict")\n'
                    f"     (SPREAD\n"
                    f'      (CASE "{pvs_traj_eq}")\n'
                    f"      ((THEN (REPLACE -1) (HIDE -1) (SKEEP)\n"
                    f"        (SPREAD (DERIV)\n"
                    f'         ((THEN (TYPEPRED "x!1") (EXPAND "{domain_type}" -1) (ASSERT))\n'
                    f'          (THEN (LEMMA "{domain_type}_dd") (INST -1 "{domain_split}")))))\n'
                    f'       (THEN (DECOMPOSE-EQUALITY 1) (HIDE 2) (TYPEPRED "x!1")\n'
                    f'        (EXPAND "{domain_type}" -1) (ASSERT)))))'
                )
            elif i == len(cases) - 1:
                # Last domain: may have piecewise splits
                if piecewise_splits:
                    # DERIVABLE branch with piecewise split
                    split_var = "x!1"
                    split_value = piecewise_splits[0]
                    proof_branches.append(
                        f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict")\n'
                        f"     (SPREAD\n"
                        f'      (CASE "{pvs_traj_eq}")\n'
                        f"      ((THEN (REPLACE -1) (DERIVABLE))\n"
                        f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "{split_var}") (EXPAND "{domain_type}" -1)\n'
                        f"        (ASSERT) (HIDE 2)\n"
                        f'        (SPREAD (CASE "{split_var}={split_value}") ((THEN (ASSERT) (GRIND)) (ASSERT)))))))'
                    )

                    # DERIV branch with piecewise split
                    split_var = "x!2"
                    proof_branches.append(
                        f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict") (SKEEP)\n'
                        f"     (SPREAD\n"
                        f'      (CASE "{pvs_traj_eq}")\n'
                        f"      ((THEN (REPLACE -1) (DERIV))\n"
                        f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "{split_var}") (EXPAND "{domain_type}")\n'
                        f'        (SPREAD (CASE "{split_var}={split_value}")\n'
                        f"         ((THEN (REPLACE -1) (ASSERT) (HIDE 2) (EVAL-EXPR 1) (ASSERT))\n"
                        f"          (ASSERT)))))))"
                    )
                else:
                    # Last domain without piecewise splits: simpler structure
                    proof_branches.append(
                        f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict")\n'
                        f"     (SPREAD\n"
                        f'      (CASE "{pvs_traj_eq}")\n'
                        f"      ((THEN (REPLACE -1) (DERIVABLE))\n"
                        f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!1") (EXPAND "{domain_type}" -1)\n'
                        f"        (ASSERT)))))"
                    )
            else:
                # Middle domains: standard CASE structure
                proof_branches.append(
                    f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict")\n'
                    f"     (SPREAD\n"
                    f'      (CASE "{pvs_traj_eq}")\n'
                    f"      ((THEN (REPLACE -1) (HIDE -1) (DERIVABLE))\n"
                    f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!1") (EXPAND "{domain_type}" -1) (FLATTEN)\n'
                    f"        (ASSERT)))))"
                )

                proof_branches.append(
                    f'(THEN (HIDE-ALL-BUT 1) (SKEEP) (EXPAND "restrict" 1) (EXPAND "g_1")\n'
                    f"     (SPREAD\n"
                    f'      (CASE "{pvs_traj_eq}")\n'
                    f'      ((THEN (REPLACE -1) (HIDE -1) (DERIV) (TYPEPRED "x!1") (EXPAND "{domain_type}" -1)\n'
                    f"        (PROPAX))\n"
                    f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!2") (EXPAND "{domain_type}" -1) (FLATTEN)\n'
                    f"        (ASSERT)))))"
                )

                proof_branches.append(
                    f'(THEN (HIDE-ALL-BUT 1) (SKEEP) (EXPAND "g_1") (EXPAND "restrict")\n'
                    f"     (SPREAD\n"
                    f'      (CASE "{pvs_traj_eq}")\n'
                    f'      ((THEN (REPLACE -1) (HIDE -1) (DERIV) (TYPEPRED "x!1") (EXPAND "{domain_type}" -1)\n'
                    f"        (PROPAX))\n"
                    f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!2") (EXPAND "{domain_type}" -1) (ASSERT)\n'
                    f"        (FLATTEN) (ASSERT)))))"
                )

    # Determine SKOLETIN argument
    skoletin_arg = "1" if num_cases >= 3 else "*"

    # Build the proof
    # Format PROPAX calls to match original: split across lines (8 per line for 3-case)
    if num_cases == 3:
        # For 3-case: 16 PROPAX split as 8 + 8
        propax_line1 = " ".join(["(PROPAX)"] * 8)
        propax_line2 = " ".join(["(PROPAX)"] * 8)
        propax_calls = f"{propax_line1}\n%|-    {propax_line2}"
    elif num_cases == 2:
        # For 2-case: 10 PROPAX, keep on one line
        propax_calls = " ".join(["(PROPAX)"] * num_propax)
    else:
        # For 4+ cases: split into groups of 8
        propax_lines = []
        remaining = num_propax
        while remaining > 0:
            count = min(8, remaining)
            propax_lines.append(" ".join(["(PROPAX)"] * count))
            remaining -= count
        propax_calls = "\n%|-    ".join(propax_lines)

    # Format each branch - don't use format_pvs_proof_line since we have explicit line breaks
    # Just add the %|- prefix to continuation lines
    formatted_branches = []
    for branch in proof_branches:
        lines = branch.split("\n")
        formatted = lines[0]  # First line doesn't need prefix
        for line in lines[1:]:
            if line.strip():
                formatted += f"\n%|-    {line}"
            else:
                formatted += "\n"
        formatted_branches.append(formatted)
    proof_branches_str = "\n%|-    ".join(formatted_branches)

    # Add trailing PROPAX calls for 3+ cases (6 trailing PROPAX for 3-case)
    if num_cases == 3:
        trailing_propax = " ".join(["(PROPAX)"] * 6)
        proof_branches_str += f"\n%|-    {trailing_propax})"
    elif num_cases > 3:
        # For 4+ cases, add some trailing PROPAX (can be adjusted based on needs)
        trailing_propax = " ".join(["(PROPAX)"] * 6)
        proof_branches_str += f"\n%|-    {trailing_propax}"

    proof = f"""%|- full_domain_soundness_lemma : PROOF
%|- (THEN (SKEEP) (SKOLETIN {skoletin_arg}) (FLATTEN) (LEMMA "full_domain_soundness_lemma_helper")
%|-  (INST -1 "x" "xo" "yo" "g_1")
%|-  (SPREAD (SPLIT -1)
%|-   ({propax_calls}
%|-    {proof_branches_str}))
%|- QED full_domain_soundness_lemma"""

    return proof


def generate_unifying_proof_from_legacy_params(
    domains: List[str],
    trajectories: List[str],
    full_traj: str,
    piecewise_split_bools: List[bool],
    domain_splits: List[str],
    piecewise_splits_map: Optional[Dict[int, List[str]]] = None,
) -> str:
    """
    Convenience wrapper that converts legacy parameter format to the new modular format.

    This function provides backward compatibility with the old function signatures.

    Args:
        domains: List of domain definitions (e.g., ["left_open(0)", "right_open(0)"])
        trajectories: List of trajectory pieces for each domain
        full_traj: Full trajectory expression
        piecewise_split_bools: Boolean list indicating which domain_splits are piecewise boundaries
        domain_splits: List of domain split points
        piecewise_splits_map: Optional dict mapping case index to piecewise splits.
            If None, will be inferred from piecewise_split_bools and domain_splits.

    Returns:
        PVS proof script as string
    """
    # Convert to case dictionaries
    # Note: In the original code, domain_splits[0] is used for the first domain
    # Other domains may not explicitly use domain_splits in the proof structure
    cases = []
    for i, (domain, trajectory) in enumerate(zip(domains, trajectories)):
        domain_type = domain.split("(")[0]
        # First domain uses domain_splits[0] if available
        domain_split = domain_splits[0] if (i == 0 and len(domain_splits) > 0) else ""

        case = {
            "domain": domain,
            "trajectory": trajectory,
            "domain_type": domain_type,
            "domain_split": domain_split,
        }
        cases.append(case)

    # Build piecewise_splits_map if not provided
    # In the original code, piecewise splits are collected globally and only used for the last domain
    if piecewise_splits_map is None:
        piecewise_splits = []
        for i in range(len(domain_splits)):
            if i < len(piecewise_split_bools) and piecewise_split_bools[i]:
                piecewise_splits.append(domain_splits[i])

        # In 3-case proof, piecewise splits are only used for the last domain (index num_cases - 1)
        if piecewise_splits and len(cases) >= 3:
            piecewise_splits_map = {len(cases) - 1: piecewise_splits}
        else:
            piecewise_splits_map = {}

    # Use modular function
    use_case_statements = len(cases) >= 3
    return generate_unifying_proof(
        cases, full_traj, piecewise_splits_map, use_case_statements
    )


# Backward compatibility functions that match old signatures
def generate_two_case_unifying_proof(
    domain_1: str,
    domain_2: str,
    trajectory_1: str,
    trajectory_2: str,
    full_traj: str,
    piecewise_split_bools: List[bool],
    domain_splits: List[str],
) -> str:
    """
    Generate proof for two cases (backward compatibility wrapper).

    This function maintains the same signature as the original in pvs_utils.py
    but uses the modular implementation.
    """
    return generate_unifying_proof_from_legacy_params(
        domains=[domain_1, domain_2],
        trajectories=[trajectory_1, trajectory_2],
        full_traj=full_traj,
        piecewise_split_bools=piecewise_split_bools,
        domain_splits=domain_splits,
    )


def generate_three_case_unifying_proof(
    domain_1: str,
    domain_2: str,
    domain_3: str,
    trajectory_1: str,
    trajectory_2: str,
    trajectory_3: str,
    full_traj: str,
    piecewise_split_bools: List[bool],
    domain_splits: List[str],
) -> str:
    """
    Generate proof for three cases (backward compatibility wrapper).

    This function maintains the same signature as the original in pvs_utils.py
    but uses the modular implementation.
    """
    return generate_unifying_proof_from_legacy_params(
        domains=[domain_1, domain_2, domain_3],
        trajectories=[trajectory_1, trajectory_2, trajectory_3],
        full_traj=full_traj,
        piecewise_split_bools=piecewise_split_bools,
        domain_splits=domain_splits,
    )


def generate_unifying_lemma_helper(
    cases: List[Dict[str, str]],
    domain_splits: List[str],
) -> str:
    """
    Generate a unifying lemma helper proof for an arbitrary number of cases.

    This is a modular replacement for generate_two_case_unifying_lemma_helper and
    generate_three_case_unifying_lemma_helper that supports any number of cases.

    Args:
        cases: List of case dictionaries, each containing:
            - lemma_name: Name of the lemma for this case (e.g., "le_lo_case_0")
            - domain: Domain definition (e.g., "left_open(0)")
            - trajectory: Trajectory function for this domain (e.g., "LAMBDA (x: real): ...")
            - case_index: Index of this case (0, 1, 2, ...) - used for EXPAND "f0", "f1", etc.
        domain_splits: List of domain split points (e.g., ["0", "4"])
            For n cases, there should be n-1 splits

    Returns:
        PVS proof script as string

    Example:
        cases = [
            {
                "lemma_name": "le_lo_case_0",
                "domain": "left_open(0)",
                "trajectory": "LAMBDA (x: real): COND x > 0 -> g(0), ELSE -> g(x) ENDCOND",
                "case_index": 0
            },
            {
                "lemma_name": "ge_ro_case_1",
                "domain": "right_open(0)",
                "trajectory": "LAMBDA (x: real): COND x >= 0 -> g(x), ELSE -> g(0) ENDCOND",
                "case_index": 1
            }
        ]
        proof = generate_unifying_lemma_helper(cases, domain_splits=["0"])
    """
    num_cases = len(cases)

    if num_cases == 0:
        return "% No cases provided"

    if num_cases == 1:
        # Single case doesn't need a helper lemma
        return "% Single case - no helper lemma needed"

    # For 3 cases, use the hardcoded version that matches the original exactly
    if num_cases == 3 and len(domain_splits) >= 2:
        return generate_three_case_unifying_lemma_helper(
            domain_split_1=domain_splits[0],
            domain_split_2=domain_splits[1],
            lemma_1=cases[0]["lemma_name"],
            lemma_2=cases[1]["lemma_name"],
            lemma_3=cases[2]["lemma_name"],
            domain_type_1=cases[0]["domain"],
            domain_type_2=cases[1]["domain"],
            domain_type_3=cases[2]["domain"],
            trajectory_function_1=cases[0]["trajectory"],
            trajectory_function_2=cases[1]["trajectory"],
            trajectory_function_3=cases[2]["trajectory"],
        )

    def generate_case_branch(case_idx: int) -> str:
        """Generate the proof branch for a single case."""
        case = cases[case_idx]
        lemma_name = case["lemma_name"]
        domain = case["domain"]
        trajectory = case["trajectory"]
        case_index = case.get("case_index", case_idx)

        # Build trajectory equality statement
        # For first case: g = trajectory_function
        # For second case: g = trajectory_function (same format as first)
        # For third case: g = trajectory_function (same format as first)
        # Note: All cases use single-line format to avoid PVS syntax errors
        traj_eq = f"restrict[real, ({domain}), real](g) = (restrict[real, ({domain}), real]({trajectory}))"

        # Determine PROPAX/ASSERT pattern based on case index
        # Pattern from original: 5 PROPAX for first case, 7 ASSERT for second, 7 ASSERT for third
        if case_idx == 0:
            propax_line = "(PROPAX) (PROPAX) (PROPAX) (PROPAX) (PROPAX)"
        elif case_idx == 1:
            # Second case: 9 ASSERT statements (from original: 7 + 2 extra)
            propax_line = "(ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT)\n%|-           (ASSERT) (ASSERT)"
        else:  # case_idx >= 2
            propax_line = (
                "(ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT) (ASSERT)"
            )

        expand_fn = f"f{case_index}"

        # Build the case proof structure
        if case_idx == 0:
            # First case structure
            # Match original formatting: CASE statement split across lines for long strings
            # The original splits long trajectory equality strings across lines
            # Add %|- prefix to all lines
            # Split traj_eq if it's long (like the original does)
            if len(traj_eq) > 80:
                # Split the trajectory equality string like the original
                # Find a good split point (after "real](g) =")
                split_point = traj_eq.find("real](g) =")
                if split_point > 0:
                    traj_eq_part1 = traj_eq[: split_point + len("real](g) =")]
                    traj_eq_part2 = traj_eq[
                        split_point + len("real](g) =") + 1 :
                    ]  # Skip space
                    return f"""%|-    (THEN (LEMMA "{lemma_name}") (INST -1 "xo" "yo" "g") (ASSERT) (EXPAND "{expand_fn}")
%|-     (ASSERT)
%|-     (SPREAD
%|-      (CASE
%|-          "{traj_eq_part1}
%|-                 {traj_eq_part2}")
%|-      ((SPREAD (SPLIT -2)
%|-        ({propax_line} (ASSERT) (ASSERT)
%|-         (THEN (INST 1 "x") (ASSERT))))
%|-       (THEN (DECOMPOSE-EQUALITY 1) (HIDE-ALL-BUT 1) (TYPEPRED "x!1")
%|-        (GRIND)))))"""
            else:
                return f"""%|-    (THEN (LEMMA "{lemma_name}") (INST -1 "xo" "yo" "g") (ASSERT) (EXPAND "{expand_fn}")
%|-     (ASSERT)
%|-     (SPREAD
%|-      (CASE
%|-          "{traj_eq}")
%|-      ((SPREAD (SPLIT -2)
%|-        ({propax_line} (ASSERT) (ASSERT)
%|-         (THEN (INST 1 "x") (ASSERT))))
%|-       (THEN (DECOMPOSE-EQUALITY 1) (HIDE-ALL-BUT 1) (TYPEPRED "x!1")
%|-        (GRIND)))))"""
        elif case_idx == num_cases - 1:
            # Last case structure
            if case_idx == 2:
                # Third case has extra ASSERT in INST line and different traj_eq format
                # For third case, traj_eq should be g = trajectory (like first case)
                # Note: case_2 ends with ))) (3 closes) not )))) (4 closes) because it's inside a list
                traj_eq_3 = f"restrict[real, ({domain}), real](g) = (restrict[real, ({domain}), real]({trajectory}))"
                return f"""%|-      (THEN (LEMMA "{lemma_name}") (INST -1 "xo" "yo" "g") (ASSERT)
%|-       (EXPAND "{expand_fn}")
%|-       (SPREAD
%|-        (CASE
%|-            "{traj_eq_3}")
%|-        ((SPREAD (SPLIT -2)
%|-          ({propax_line}
%|-           (THEN (ASSERT) (INST 1 "x") (ASSERT))))
%|-         (THEN (HIDE-ALL-BUT 1) (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!1")
%|-          (GRIND)))))"""
            else:
                return f"""%|-      (THEN (LEMMA "{lemma_name}") (INST -1 "xo" "yo" "g") (ASSERT)
%|-       (EXPAND "{expand_fn}")
%|-       (SPREAD
%|-        (CASE
%|-            "{traj_eq}")
%|-        ((SPREAD (SPLIT -2)
%|-          ({propax_line}
%|-           (THEN (INST 1 "x") (ASSERT))))
%|-         (THEN (HIDE-ALL-BUT 1) (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!1")
%|-          (GRIND))))))"""
        else:
            # Middle cases (for 4+ cases) - same structure as case 2
            return f"""%|-      (THEN (LEMMA "{lemma_name}") (INST -1 "xo" "yo" "g") (ASSERT)
%|-       (EXPAND "{expand_fn}")
%|-       (SPREAD
%|-        (CASE
%|-            "{traj_eq}")
%|-        ((SPREAD (SPLIT -2)
%|-          ({propax_line}
%|-           (THEN (INST 1 "x") (ASSERT))))
%|-         (THEN (ASSERT) (DECOMPOSE-EQUALITY 1) (HIDE-ALL-BUT 1) (TYPEPRED "x!1")
%|-          (GRIND))
%|-         (THEN (ASSERT) (HIDE-ALL-BUT 1) (GRIND)))))"""

    # Systematic approach: strip all prefixes, build structure, add prefixes back once
    # Helper: strip %|- prefix from all lines recursively
    def strip_all_prefixes(text):
        """Remove %|- prefix from all lines, keeping content and relative indentation."""
        lines = text.split("\n")
        result = []
        for line in lines:
            # Remove all occurrences of %|- prefix (handles duplicates)
            cleaned = line
            while "%|-" in cleaned and cleaned.strip().startswith("%|-"):
                idx = cleaned.find("%|-")
                if idx != -1:
                    cleaned = cleaned[idx + 3 :]
                else:
                    break
            result.append(cleaned)
        return "\n".join(result)

    # Add %|- prefix to all non-empty lines, preserving indentation
    def add_prefixes(text):
        """Add %|- prefix to all non-empty lines, preserving existing indentation."""
        lines = text.split("\n")
        result = []
        for line in lines:
            if line.strip():
                # Count leading spaces
                leading_spaces = len(line) - len(line.lstrip())
                # Preserve indentation: %|- followed by the same number of spaces
                if leading_spaces > 0:
                    result.append(f"%|-{' ' * leading_spaces}{line.lstrip()}")
                else:
                    result.append(f"%|-{line}")
            else:
                result.append("")
        return "\n".join(result)

    # Generate all case branches and strip prefixes
    case_branches_raw = [
        strip_all_prefixes(generate_case_branch(i)) for i in range(num_cases)
    ]

    # Build nested structure systematically for all case counts
    if num_cases == 2:
        # Two cases: simple structure
        nested_structure = f"""   ({case_branches_raw[0]}
    {case_branches_raw[1]})"""
    elif num_cases == 3:
        # Three cases: one level of nesting - match original formatting exactly
        # Structure: (SPREAD (CASE "x <= 0")
        #   ((THEN case_0 ...)))  [case_0 ends with ))) closing its own structure]
        #   (SPREAD (CASE "x <= 4")
        #     ((THEN case_1 ...)))  [case_1 ends with ))) closing its own structure]
        #     (THEN case_2 ...)))  [case_2 ends with ))) closing: GRIND's THEN, CASE list, SPREAD list]
        #     ))  [close inner SPREAD's list and SPREAD]
        #   ))  [close outer SPREAD's list and SPREAD]
        # ))  [close outer CASE and THEN]
        # Total: case_2's ))) (3) + inner )) (2) + outer )) (2) + final )) (2) = 9 closes after GRIND
        # Strip leading spaces from all case branches since we're providing the indentation
        # Also ensure they start with (THEN, not ((THEN
        case_0_stripped = case_branches_raw[0].lstrip()
        if case_0_stripped.startswith("(("):
            case_0_stripped = case_0_stripped[1:]  # Remove extra opening paren
        case_1_stripped = case_branches_raw[1].lstrip()
        if case_1_stripped.startswith("(("):
            case_1_stripped = case_1_stripped[1:]  # Remove extra opening paren
        case_2_stripped = case_branches_raw[2].lstrip()
        # Build structure matching original exactly
        # The structure is:
        #   (SPREAD (CASE "x <= split_1")
        #    ((THEN case_0 ...)))
        #    (SPREAD (CASE "x <= split_2")
        #     ((THEN case_1 ...)))
        #     (THEN case_2 ...))))))
        #
        # case_0_stripped already starts with (THEN, so we only add one ( before it
        # case_1_stripped already starts with (THEN, so we only add one ( before it
        nested_structure = f"""   ({case_0_stripped}
    (SPREAD (CASE "x <= {domain_splits[1]}")
     ({case_1_stripped}
      {case_2_stripped}))))))"""
    else:
        # Four or more cases: build nested structure iteratively
        nested_cases_raw = case_branches_raw[1]
        for i in range(2, num_cases):
            if i - 1 < len(domain_splits):
                nested_cases_raw = f"""    (SPREAD (CASE "x <= {domain_splits[i-1]}")
     (
{nested_cases_raw}
      {case_branches_raw[i]}))"""
            else:
                nested_cases_raw = f"""{nested_cases_raw}
      {case_branches_raw[i]}"""
        nested_structure = f"""   ({case_branches_raw[0]}
    {nested_cases_raw})"""

    # Add prefixes once to the entire structure
    nested_prefixed = add_prefixes(nested_structure)

    # Build final proof - nested_prefixed already has all prefixes
    # For 3-case, nested_structure already includes all 9 closing parens after case_2's GRIND
    # For other cases, we need to add closing parens here
    if num_cases == 3:
        # All closing parens (9 total) are already in nested_structure
        closing_parens = ""
    elif num_cases == 2:
        # For 2-case, need to close the structure
        closing_parens = "))"
    else:
        # For 4+ cases, add closing parens
        closing_parens = "))"
    proof = f"""%|- full_domain_soundness_lemma_helper : PROOF
%|- (THEN (SKEEP) (SKOLETIN*) (FLATTEN)
%|-  (SPREAD (CASE "x <= {domain_splits[0]}")
{nested_prefixed}{closing_parens}
%|- QED full_domain_soundness_lemma_helper"""

    return proof


# Backward compatibility functions for helper proofs
def generate_two_case_unifying_lemma_helper(
    domain_split: str,
    lemma_1: str,
    lemma_2: str,
    domain_type_1: str,
    domain_type_2: str,
    trajectory_function_1: str,
    trajectory_function_2: str,
) -> str:
    """
    Generate helper proof for two cases (backward compatibility wrapper).

    This function maintains the same signature as the original in pvs_utils.py
    but uses the modular implementation.
    """
    cases = [
        {
            "lemma_name": lemma_1,
            "domain": domain_type_1,
            "trajectory": trajectory_function_1,
            "case_index": 0,
        },
        {
            "lemma_name": lemma_2,
            "domain": domain_type_2,
            "trajectory": trajectory_function_2,
            "case_index": 1,
        },
    ]
    return generate_unifying_lemma_helper(cases, domain_splits=[domain_split])


def generate_three_case_unifying_lemma_helper(
    domain_split_1: str,
    domain_split_2: str,
    lemma_1: str,
    lemma_2: str,
    lemma_3: str,
    domain_type_1: str,
    domain_type_2: str,
    domain_type_3: str,
    trajectory_function_1: str,
    trajectory_function_2: str,
    trajectory_function_3: str,
) -> str:
    """
    Generate helper proof for three cases (backward compatibility wrapper).

    This function maintains the same signature as the original in pvs_utils.py
    and matches the original formatting exactly.
    """
    # For 3-case, use the exact original format to ensure compatibility
    # Split long trajectory equality strings like the original does
    traj_eq_1 = f"restrict[real, ({domain_type_1}), real](g) = (restrict[real, ({domain_type_1}), real]"
    traj_eq_1_cont = f"                 ({trajectory_function_1}))"

    return f"""%|- full_domain_soundness_lemma_helper : PROOF
%|- (THEN (SKEEP) (SKOLETIN*) (FLATTEN)
%|-  (SPREAD (CASE "x <= {domain_split_1}")
%|-   ((THEN (LEMMA "{lemma_1}") (INST -1 "xo" "yo" "g") (ASSERT) (EXPAND "f0")
%|-     (ASSERT)
%|-     (SPREAD
%|-      (CASE
%|-          "{traj_eq_1}
%|-                 {traj_eq_1_cont}")
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
%|- QED full_domain_soundness_lemma_helper
"""
