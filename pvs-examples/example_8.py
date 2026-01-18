#!/usr/bin/env python
# coding: utf-8
"""
4-case trajectory example to test proof_builder extensibility.

f(x) =
  x²               for x ≤ -10     (parabola)
  -20x - 100       for -10 < x ≤ 0 (linear)
  2x² - 20x - 100  for 0 < x ≤ 10  (parabola)
  20x - 300        for x > 10      (linear)

This function is C¹ continuous.
Domains: left_open(-10), ci(-10, 0), ci(0, 10), right_open(10)

Each interval has width ≥ 10, which is larger than the polygon width of 4,
ensuring the "fully inside" case is not spurious.

C¹ continuity verification:
  At x = -10: f₁(-10) = 100, f₂(-10) = 100; f₁'(-10) = -20, f₂'(-10) = -20 ✓
  At x = 0:   f₂(0) = -100, f₃(0) = -100; f₂'(0) = -20, f₃'(0) = -20 ✓
  At x = 10:  f₃(10) = -100, f₄(10) = -100; f₃'(10) = 20, f₄'(10) = 20 ✓
"""

import argparse
from sympy import symbols, Point, Polygon, Interval, oo, Piecewise, And, plot
from pvs_utils import (
    generate_complete_proof_package,
    log_proof_to_file,
)


def main():
    parser = argparse.ArgumentParser(description="4-case trajectory example")
    parser.add_argument(
        "output_filename",
        nargs="?",
        default="example_8.pvs",
        help="Output filename for the proof (default: example_8.pvs)",
    )

    args = parser.parse_args()

    x = symbols("x")

    # Create a polygon (rectangle)
    w = 1.0
    rect_points = [
        Point(val) for val in [[2 * w, -w], [2 * w, w], [-2 * w, w], [-2 * w, -w]]
    ]
    polygon = Polygon(*rect_points)

    # 4-case piecewise trajectory with wider intervals (each > 4)
    # C¹ continuity at each breakpoint (-10, 0, 10)
    trajectory_expr = Piecewise(
        (x**2, x <= -10),  # parabola: f(-10)=100, f'(-10)=-20
        (-20 * x - 100, And(x > -10, x <= 0)),  # linear: f(0)=-100, f'(0)=-20
        (2 * x**2 - 20 * x - 100, And(x > 0, x <= 10)),  # parabola: f(10)=-100, f'(10)=20
        (20 * x - 300, x > 10),  # linear: matches at x=10
    )
    # plot(trajectory_expr)
    # trajectory_expr = Piecewise((x**2 + 2 * x + 1, x <= 3), (8 * x - 8, x > 3))

    domain = Interval(-oo, oo)

    package = generate_complete_proof_package(trajectory_expr, polygon, domain)

    log_proof_to_file(package, args.output_filename)

    print(f"Proof generated successfully and saved to: {args.output_filename}")


if __name__ == "__main__":
    main()
