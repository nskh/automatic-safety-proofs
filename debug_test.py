#!/usr/bin/env python
# coding: utf-8
"""
Debug test for automatic proof generation.
"""

from sympy import *
from pvs_utils import generate_unbounded_proof_calls, generate_proof_scripts_from_calls

# Define symbols
x, y = symbols("x y")

# Create a square polygon with width 2 (half-width 1)
w = 1.0
square_points = [Point(val) for val in [[w, -w], [w, w], [-w, w], [-w, -w]]]
square = Polygon(*square_points)

# Use y = x^2 trajectory which should have transition point at (0,0)
trajectory = x**2  # This should have transition point at x=0 where derivative is 0

# Domain (infinite)
domain = Interval(-oo, oo)

print("Debug test with y = x^2 trajectory...")
print(f"Trajectory: {trajectory}")
print(f"Domain: {domain}")

# Add debug prints to understand the transition points
print(f"\nDomain bounds: inf={domain.inf}, sup={domain.sup}")

try:
    # Generate proof calls
    proof_calls = generate_unbounded_proof_calls(trajectory, square, domain)

    print(f"\nGenerated {len(proof_calls)} proof calls:")
    for i, call in enumerate(proof_calls):
        print(f"\nProof call {i+1}:")
        print(f"  Case label: {call['case_label']}")
        print(f"  Deriv lemma: {call['deriv_lemma']}")
        print(f"  Max right: {call['max_right']}")
        print(f"  Min left: {call['min_left']}")
        print(f"  Domain definition: {call['domain_definition']}")
        print(f"  Interval: {call['interval']}")
        print(f"  Active corner: {call['active_corner']}")

    # Generate proof scripts
    print(f"\n{'='*60}")
    print("GENERATING PROOF SCRIPTS:")
    print(f"{'='*60}")

    proof_scripts = generate_proof_scripts_from_calls(proof_calls)
    print(f"\nGenerated {len(proof_scripts)} proof scripts:")

    for i, script in enumerate(proof_scripts):
        print(f"\n{'='*60}")
        print(f"PROOF SCRIPT {i+1}:")
        print(f"{'='*60}")
        print(script)
        print(f"{'='*60}")

except Exception as e:
    print(f"Error: {e}")
    import traceback

    traceback.print_exc()
