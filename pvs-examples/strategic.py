#!/usr/bin/env python
# coding: utf-8
"""
Strategic example with three-piece linear piecewise trajectory.
"""

import argparse
from sympy import symbols, Point, Polygon, RegularPolygon, Interval, oo, Piecewise
from pvs_utils import (
    generate_complete_proof_package,
    log_proof_to_file,
)


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Example for automatic proof generation"
    )
    parser.add_argument(
        "output_filename",
        nargs="?",
        default="strategic.pvs",
        help="Output filename for the proof (default: strategic.pvs)",
    )

    args = parser.parse_args()

    # Define symbols
    x = symbols("x")

    # Create a 2-1 rectangle (standard from other examples)
    w = 1.0
    rect_points = [
        Point(val) for val in [[2 * w, -w], [2 * w, w], [-2 * w, w], [-2 * w, -w]]
    ]
    polygon = Polygon(*rect_points)

    # Trajectory function: three-piece linear piecewise
    # Intervals are [start, end) - each piece includes its starting boundary
    trajectory_expr = Piecewise(
        (73 * x - 3011, x < 4),
        (81 * x - 3043, (x >= 4) & (x < 21)),
        (71 * x - 2833, x >= 21),
    )

    # Domain
    domain = Interval(-oo, oo)

    package = generate_complete_proof_package(trajectory_expr, polygon, domain)

    log_proof_to_file(package, args.output_filename)

    print(f"Proof generated successfully and saved to: {args.output_filename}")


if __name__ == "__main__":
    main()
