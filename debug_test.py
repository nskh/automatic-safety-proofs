#!/usr/bin/env python
# coding: utf-8
"""
Debug test for automatic proof generation with lemma generation.
"""

from sympy import symbols, Point, Polygon, RegularPolygon, Interval, oo, Piecewise
from pvs_utils import (
    generate_complete_proof_package,
    print_prooflite,
    log_proof_to_file,
)

# Define symbols
x = symbols("x")

# Create a rectangle polygon with width 4 (half-width 2)
w = 1.0
rect_points = [
    Point(val) for val in [[2 * w, -w], [2 * w, w], [-2 * w, w], [-2 * w, -w]]
]
polygon = Polygon(*rect_points)
# diamond_points = [Point(val) for val in [[0, 1], [1, 0], [0, -1], [-1, 0]]]
# polygon = Polygon(*diamond_points)
# polygon = RegularPolygon((0, 0), 1, n=6)  # hexagon

# trajectory_expr = x**2
trajectory_expr = Piecewise((0, x <= 4), (8 * x - 16, x > 4))
# trajectory_expr = Piecewise((0, (x < 3.0)), (1000671.88367391 - 44141.8652246365*x, (x < 26.0)), (292190.671499991 - 16892.587833333*x, (x <= 27.0)), (0, x>27.0))
# trajectory_expr = Piecewise((-2.385292, (x < 3.0)), (0.121268860507243*x - 2.74909858152173, (x >= 3.0) & (x < 26.0)), (0.0464082083333324*x - 0.802721624999975, (x >= 26.0) & (x <= 27.0)), (0.0464082083333324*27.0 - 0.802721624999975, (x>27.0)))
# trajectory_expr = Piecewise((-2.3853, (x < 3.0) & (x >= 0)), (0.1213*x - 2.7491, (x >= 3.0) & (x < 26.0)), (0.0464*x - 0.8027, (x >= 26.0) & (x <= 27.0)))
# trajectory_expr = Piecewise((-2.385292, x<= 3.0), (0.121268860507243*x - 2.74909858152173, (x > 3.0) & (x < 26.0)), (0.0464082083333324*x - 0.802721624999975, (x >= 26.0) & (x <= 27.0)), (0.0464082083333324*27.0 - 0.802721624999975, x>27.0))
# trajectory_expr = Piecewise((0.121268860507243*x - 2.74909858152173, (x < 26.0)), (0.0464082083333324*x - 0.802721624999975, (x >= 26.0)))
# trajectory_expr = Piecewise((-2.385292, (x < 3.0)), (0.121268860507243*x - 2.74909858152173, (x >= 3.0)))
# trajectory_expr = Piecewise((- 1505533.5482712, (x <0)),(39741.5279834603*x - 1505533.5482712, (x < 1.0) & (x >= 0)), (32975.9468185201*x - 1498767.96710626, (x >= 1.0) & (x < 3.0)), (41270.9129515472*x - 1523652.86550535, (x >= 3.0) & (x < 4.0)), (41183.7884000142*x - 1523304.36729921, (x >= 4.0) & (x < 5.0)), (41096.458085676*x - 1522867.71572752, (x >= 5.0) & (x < 6.0)), (41008.9223106166*x - 1522342.50107717, (x >= 6.0) & (x < 7.0)), (40921.1813773413*x - 1521728.31454424, (x >= 7.0) & (x < 8.0)), (40833.2355890456*x - 1521024.74823787, (x >= 8.0) & (x < 9.0)), (40745.085249475*x - 1520231.39518174, (x >= 9.0) & (x < 10.0)), (40656.7306629571*x - 1519347.84931656, (x >= 10.0) & (x < 11.0)), (40568.1721343962*x - 1518373.70550239, (x >= 11.0) & (x < 12.0)), (40479.4099692759*x - 1517308.55952095, (x >= 12.0) & (x < 13.0)), (40390.4444736562*x - 1516152.00807789, (x >= 13.0) & (x < 14.0)), (40301.2759541732*x - 1514903.64880513, (x >= 14.0) & (x < 15.0)), (40211.9047180382*x - 1513563.0802631, (x >= 15.0) & (x < 16.0)), (40122.3310730384*x - 1512129.90194311, (x >= 16.0) & (x < 17.0)), (40032.5553275337*x - 1510603.71426953, (x >= 17.0) & (x < 18.0)), (39942.5777904574*x - 1508984.11860215, (x >= 18.0) & (x < 19.0)), (39852.3987713151*x - 1507270.71723845, (x >= 19.0) & (x < 20.0)), (39762.0185801843*x - 1505463.11341583, (x >= 20.0) & (x < 21.0)), (39671.4375277126*x - 1503560.91131393, (x >= 21.0) & (x < 22.0)), (39205.8585682418*x - 1493318.17420557, (x >= 22.0) & (x < 23.0)), (38454.5226300755*x - 1476037.44762774, (x >= 23.0) & (x < 24.0)), (38511.8183428739*x - 1477412.54473491, (x >= 24.0) & (x < 25.0)), (38569.0377449836*x - 1478843.02978765, (x >= 25.0) & (x < 26.0)), (19139.4993661304*x - 973675.031937467, (x >= 26.0) & (x <= 27.0)), (19139.4993661304*27.0 - 973675.031937467, (x>27.0)))
#Piecewise((7000.0, (x < 1.0) & (x >= 0)), (8000.0 - 1000.0*x, (x >= 1.0) & (x < 2.0)), (6000.0, (x >= 2.0) & (x < 3.0)), (9000.0 - 1000.0*x, (x >= 3.0) & (x < 4.0)), (5000.0, (x >= 4.0) & (x < 20.0)), (25000.0 - 1000.0*x, (x >= 20.0) & (x < 21.0)), (4000.0, (x >= 21.0) & (x < 26.0)), (30000.0 - 1000.0*x, (x >= 26.0) & (x <= 27.0)))
# Domain (infinite)
domain = Interval(-oo, oo)
# domain = Interval(3, oo)
# domain = Interval(-oo, 10)
# domain = Interval(0, 27)

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
        trajectory_expr, polygon, domain, "debug_lemma"
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
        print(f"  Case labels: {call['case_labels']}")
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

    print(print_prooflite(package))

    log_proof_to_file(package, "debug_proof.pvs")

    # print(f"Debug: Number of lemmas: {len(package['lemmas'])}")
    # print(f"Debug: Number of proof scripts: {len(package['proof_scripts'])}")
    # print(f"Debug: Lemmas type: {type(package['lemmas'])}")
    # print(
    #     f"Debug: First lemma: {package['lemmas'][0] if package['lemmas'] else 'None'}"
    # )

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
