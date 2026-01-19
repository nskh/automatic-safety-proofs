#!/usr/bin/env python
# coding: utf-8
"""
4-case trajectory example similar to example_8 with different parameters.

f(x) =
  x²               for x ≤ -8      (parabola)
  -16x - 64        for -8 < x ≤ 0  (linear)
  2x² - 16x - 64   for 0 < x ≤ 8   (parabola with vertex at x=4)
  16x - 192        for x > 8       (linear)

This function is C¹ continuous.
Domains: left_open(-8), ci(-8, 0), ci(0, 4), ci(4, 8), right_open(8)
(5 domains due to notch at x=4 where derivative = 0)

C¹ continuity verification:
  At x = -8: f₁(-8) = 64, f₂(-8) = 64; f₁'(-8) = -16, f₂'(-8) = -16 ✓
  At x = 0:  f₂(0) = -64, f₃(0) = -64; f₂'(0) = -16, f₃'(0) = -16 ✓
  At x = 8:  f₃(8) = -64, f₄(8) = -64; f₃'(8) = 16, f₄'(8) = 16 ✓
"""

import argparse
from sympy import symbols, Point, Polygon, Interval, oo, Piecewise, And
from pvs_utils import (
    generate_complete_proof_package,
    log_proof_to_file,
)


def main():
    parser = argparse.ArgumentParser(description="4-case trajectory example")
    parser.add_argument(
        "output_filename",
        nargs="?",
        default="example_9.pvs",
        help="Output filename for the proof (default: example_9.pvs)",
    )

    args = parser.parse_args()

    x = symbols("x")

    # Create a polygon (rectangle)
    w = 1.0
    rect_points = [
        Point(val) for val in [[2 * w, -w], [2 * w, w], [-2 * w, w], [-2 * w, -w]]
    ]
    polygon = Polygon(*rect_points)

    # 4-case piecewise trajectory with C¹ continuity at each breakpoint (-8, 0, 8)
    # Notch at x=4 where f₃'(4) = 0
    trajectory_expr = Piecewise(
        (x**2, x <= -8),                              # parabola: f(-8)=64, f'(-8)=-16
        (-16 * x - 64, And(x > -8, x <= 0)),          # linear: f(0)=-64, f'(0)=-16
        (2 * x**2 - 16 * x - 64, And(x > 0, x <= 8)), # parabola: vertex at x=4, f(8)=-64, f'(8)=16
        (16 * x - 192, x > 8),                        # linear: matches at x=8
    )

    domain = Interval(-oo, oo)

    package = generate_complete_proof_package(trajectory_expr, polygon, domain)

    log_proof_to_file(package, args.output_filename)

    print(f"Proof generated successfully and saved to: {args.output_filename}")


if __name__ == "__main__":
    main()
