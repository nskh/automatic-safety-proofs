#!/usr/bin/env python
# coding: utf-8
"""
Debug test for automatic proof generation with lemma generation.
"""

from sympy import *
from pvs_utils import (
    generate_unbounded_proof_calls,
    generate_proof_scripts_from_calls,
    generate_complete_proof_package,
    generate_lemmas_from_calls,
)

# Define symbols
x, y = symbols("x y")

# Create a square polygon with width 2 (half-width 1)
w = 1.0
square_points = [Point(val) for val in [[w, -w], [w, w], [-w, w], [-w, -w]]]
square = Polygon(*square_points)

# Use y = x^2 trajectory which should have transition point at (0,0)
trajectory_expr = x**2  # This should have transition point at x=0 where derivative is 0

# Domain (infinite)
domain = Interval(-oo, oo)

print("Debug test with y = x^2 trajectory and lemma generation...")
print(f"Trajectory: {trajectory_expr}")
print(f"Domain: {domain}")

# Add debug prints to understand the transition points
print(f"\nDomain bounds: inf={domain.inf}, sup={domain.sup}")

try:
    # Generate complete proof package including lemmas
    print(f"\n{'='*60}")
    print("GENERATING COMPLETE PROOF PACKAGE:")
    print(f"{'='*60}")

    package = generate_complete_proof_package(
        trajectory_expr, square, domain, "debug_lemma"
    )

    print(f"Package keys: {list(package.keys())}")
    print(f"Number of proof calls: {len(package['proof_calls'])}")
    print(f"Number of lemmas: {len(package['lemmas'])}")
    print(f"Number of proof scripts: {len(package['proof_scripts'])}")

    # Display proof calls with lemma information
    print(f"\n{'='*60}")
    print("PROOF CALLS WITH LEMMA INFORMATION:")
    print(f"{'='*60}")

    for i, call in enumerate(package["proof_calls"]):
        print(f"\nProof call {i+1}:")
        print(f"  Case label: {call['case_label']}")
        print(f"  Deriv lemma: {call['deriv_lemma']}")
        print(f"  Max right: {call['max_right']}")
        print(f"  Min left: {call['min_left']}")
        print(f"  Domain definition: {call['domain_definition']}")
        print(f"  Interval: {call['interval']}")
        print(f"  Active corner: {call['active_corner']}")
        print(f"  Lemma name: {call.get('lemma_name', 'N/A')}")
        print(f"  Deriv clause: {call.get('deriv_clause1', 'N/A')}")
        print(f"  Domain start: {call.get('domain_start', 'N/A')}")
        print(f"  Domain end: {call.get('domain_end', 'N/A')}")

    # Display lemmas and proof scripts together
    print(f"\n{'='*60}")
    print("LEMMAS AND PROOF SCRIPTS TOGETHER:")
    print(f"{'='*60}")

    print(f"Debug: Number of lemmas: {len(package['lemmas'])}")
    print(f"Debug: Number of proof scripts: {len(package['proof_scripts'])}")
    print(f"Debug: Lemmas type: {type(package['lemmas'])}")
    print(
        f"Debug: First lemma: {package['lemmas'][0] if package['lemmas'] else 'None'}"
    )

    for i in range(len(package["lemmas"])):
        print(f"\n{'='*60}")
        print(f"CASE {i+1}:")
        print(f"{'='*60}")

        # Print the lemma
        print(f"\nLEMMA {i+1}:")
        print(f"{'='*40}")
        print(package["lemmas"][i])
        print(f"{'='*40}")

        # Print the corresponding proof script
        print(f"\nPROOF SCRIPT {i+1}:")
        print(f"{'='*40}")
        print(package["proof_scripts"][i])
        print(f"{'='*40}")
except Exception as e:
    print(f"Error: {e}")
    import traceback

    traceback.print_exc()
