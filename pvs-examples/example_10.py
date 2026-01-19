#!/usr/bin/env python
# coding: utf-8
"""
3-case trajectory with 4 domains: linear-parabola-linear (V-shape with curved bottom).

This is substantively different from examples 8 and 9 which are parabola-first.

f(x) =
  -6x - 72    for x ≤ -12    (linear descending, f(-12)=0, f'=-6)
  x²/4 - 36   for -12 < x ≤ 12 (parabola with notch at x=0, f(-12)=0, f(0)=-36, f(12)=0)
  6x - 72     for x > 12      (linear ascending, f(12)=0, f'=6)

This function is C¹ continuous at all breakpoints.
Domains: left_open(-12), ci(-12, 0), ci(0, 12), right_open(12)
(4 domains from 3 pieces due to notch at x=0 where derivative = 0)

C¹ continuity verification:
  At x = -12: f₁(-12) = 72-72 = 0, f₂(-12) = 36-36 = 0; f₁'(-12) = -6, f₂'(-12) = -6 ✓
  At x = 12:  f₂(12) = 36-36 = 0, f₃(12) = 72-72 = 0; f₂'(12) = 6, f₃'(12) = 6 ✓

Notch at x = 0: f₂'(x) = x/2 = 0 at x = 0
"""

import argparse
from sympy import symbols, Point, Polygon, Interval, oo, Piecewise, And
from pvs_utils import (
    generate_complete_proof_package,
    log_proof_to_file,
)


def main():
    parser = argparse.ArgumentParser(description="3-case linear-parabola-linear trajectory with 4 domains")
    parser.add_argument(
        "output_filename",
        nargs="?",
        default="example_10.pvs",
        help="Output filename for the proof (default: example_10.pvs)",
    )

    args = parser.parse_args()

    x = symbols("x")

    # Create a polygon (rectangle)
    w = 1.0
    rect_points = [
        Point(val) for val in [[2 * w, -w], [2 * w, w], [-2 * w, w], [-2 * w, -w]]
    ]
    polygon = Polygon(*rect_points)

    # 3-case piecewise: linear, parabola (with notch at x=0), linear
    # All intervals width 12 > polygon width 4
    trajectory_expr = Piecewise(
        (-6 * x - 72, x <= -12),                          # linear descending
        (x**2 / 4 - 36, And(x > -12, x <= 12)),           # parabola with vertex at x=0
        (6 * x - 72, x > 12),                             # linear ascending
    )

    domain = Interval(-oo, oo)

    package = generate_complete_proof_package(trajectory_expr, polygon, domain)

    log_proof_to_file(package, args.output_filename)

    print(f"Proof generated successfully and saved to: {args.output_filename}")


if __name__ == "__main__":
    main()
