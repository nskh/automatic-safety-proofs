"""
Modular implementation of unifying proof generation.

This module provides a flexible interface for generating PVS proof scripts
that can handle an arbitrary number of cases, replacing the fixed-case
implementations in pvs_utils.py.
"""

from typing import List, Dict, Optional, Union


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
            f'     (SPREAD (DERIV)',
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
        
        for i, (domain_type, pvs_traj_eq) in enumerate(zip(domain_types, pvs_traj_equalities)):
            domain_split = cases[i].get("domain_split", "")
            
            # Generate branches for this domain
            if i == 0:
                # First domain: DERIVABLE branch with CASE
                proof_branches.append(
                    f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict")'
                    f'     (SPREAD (DERIVABLE)'
                    f'      ((SPREAD'
                    f'        (CASE "{pvs_traj_eq}")'
                    f'        ((THEN (REPLACE -1) (HIDE -1) (DERIVABLE) (LEMMA "{domain_type}_dd")'
                    f'          (INST -1 "{domain_split}"))'
                    f'         (THEN (DECOMPOSE-EQUALITY 1) (HIDE 2) (TYPEPRED "x!1")'
                    f'          (EXPAND "{domain_type}" -1) (ASSERT))))'
                    f'       (THEN (LEMMA "{domain_type}_dd") (INST -1 "{domain_split}")))))'
                )
                
                # First domain: DERIV branch with CASE
                proof_branches.append(
                    f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict")'
                    f'     (SPREAD'
                    f'      (CASE "{pvs_traj_eq}")'
                    f'      ((THEN (REPLACE -1) (HIDE -1) (SKEEP)'
                    f'        (SPREAD (DERIV)'
                    f'         ((THEN (TYPEPRED "x!1") (EXPAND "{domain_type}" -1) (ASSERT))'
                    f'          (THEN (LEMMA "{domain_type}_dd") (INST -1 "{domain_split}")))))'
                    f'       (THEN (DECOMPOSE-EQUALITY 1) (HIDE 2) (TYPEPRED "x!1")'
                    f'        (EXPAND "{domain_type}" -1) (ASSERT)))))'
                )
            elif i == len(cases) - 1:
                # Last domain: may have piecewise splits
                if piecewise_splits:
                    # DERIVABLE branch with piecewise split
                    split_var = "x!1"
                    split_value = piecewise_splits[0]
                    proof_branches.append(
                        f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict")'
                        f'     (SPREAD'
                        f'      (CASE "{pvs_traj_eq}")'
                        f'      ((THEN (REPLACE -1) (DERIVABLE))'
                        f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "{split_var}") (EXPAND "{domain_type}" -1)'
                        f'        (ASSERT) (HIDE 2)'
                        f'        (SPREAD (CASE "{split_var}={split_value}") ((THEN (ASSERT) (GRIND)) (ASSERT)))))))'
                    )
                    
                    # DERIV branch with piecewise split
                    split_var = "x!2"
                    proof_branches.append(
                        f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict") (SKEEP)'
                        f'     (SPREAD'
                        f'      (CASE "{pvs_traj_eq}")'
                        f'      ((THEN (REPLACE -1) (DERIV))'
                        f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "{split_var}") (EXPAND "{domain_type}")'
                        f'        (SPREAD (CASE "{split_var}={split_value}")'
                        f'         ((THEN (REPLACE -1) (ASSERT) (HIDE 2) (EVAL-EXPR 1) (ASSERT))'
                        f'          (ASSERT)))))))'
                    )
                else:
                    # Last domain without piecewise splits: simpler structure
                    proof_branches.append(
                        f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict")'
                        f'     (SPREAD'
                        f'      (CASE "{pvs_traj_eq}")'
                        f'      ((THEN (REPLACE -1) (DERIVABLE))'
                        f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!1") (EXPAND "{domain_type}" -1)'
                        f'        (ASSERT)))))'
                    )
            else:
                # Middle domains: standard CASE structure
                proof_branches.append(
                    f'(THEN (HIDE-ALL-BUT 1) (EXPAND "g_1") (EXPAND "restrict")'
                    f'     (SPREAD'
                    f'      (CASE "{pvs_traj_eq}")'
                    f'      ((THEN (REPLACE -1) (HIDE -1) (DERIVABLE))'
                    f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!1") (EXPAND "{domain_type}" -1) (FLATTEN)'
                    f'        (ASSERT)))))'
                )
                
                proof_branches.append(
                    f'(THEN (HIDE-ALL-BUT 1) (SKEEP) (EXPAND "restrict" 1) (EXPAND "g_1")'
                    f'     (SPREAD'
                    f'      (CASE "{pvs_traj_eq}")'
                    f'      ((THEN (REPLACE -1) (HIDE -1) (DERIV) (TYPEPRED "x!1") (EXPAND "{domain_type}" -1)'
                    f'        (PROPAX))'
                    f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!2") (EXPAND "{domain_type}" -1) (FLATTEN)'
                    f'        (ASSERT)))))'
                )
                
                proof_branches.append(
                    f'(THEN (HIDE-ALL-BUT 1) (SKEEP) (EXPAND "g_1") (EXPAND "restrict")'
                    f'     (SPREAD'
                    f'      (CASE "{pvs_traj_eq}")'
                    f'      ((THEN (REPLACE -1) (HIDE -1) (DERIV) (TYPEPRED "x!1") (EXPAND "{domain_type}" -1)'
                    f'        (PROPAX))'
                    f'       (THEN (DECOMPOSE-EQUALITY 1) (TYPEPRED "x!2") (EXPAND "{domain_type}" -1) (ASSERT)'
                    f'        (FLATTEN) (ASSERT)))))'
                )
    
    # Determine SKOLETIN argument
    skoletin_arg = "1" if num_cases >= 3 else "*"
    
    # Build the proof
    propax_calls = " ".join(["(PROPAX)"] * num_propax)
    
    proof_branches_str = "\n%|-    ".join(proof_branches)
    
    # Add trailing PROPAX calls for 3+ cases (6 trailing PROPAX for 3-case)
    if num_cases == 3:
        trailing_propax = " ".join(["(PROPAX)"] * 6)
        proof_branches_str += f"\n%|-    {trailing_propax}"
    elif num_cases > 3:
        # For 4+ cases, add some trailing PROPAX (can be adjusted based on needs)
        trailing_propax = " ".join(["(PROPAX)"] * 6)
        proof_branches_str += f"\n%|-    {trailing_propax}"
    
    proof = f"""%|- full_domain_soundness_lemma : PROOF
%|- (THEN (SKEEP) (SKOLETIN {skoletin_arg}) (FLATTEN) (LEMMA "full_domain_soundness_lemma_helper")
%|-  (INST -1 "x" "xo" "yo" "g_1")
%|-  (SPREAD (SPLIT -1)
%|-   ({propax_calls}
%|-    {proof_branches_str})))
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
    return generate_unifying_proof(cases, full_traj, piecewise_splits_map, use_case_statements)


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

