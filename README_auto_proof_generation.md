# Automatic Proof Generation

This module provides functionality to automatically generate calls to `unbounded_one_side_proof_script` based on trajectory, polygon, and domain analysis.

## Overview

The automatic proof generation system analyzes:

1. **Trajectory**: Function of x (y=f(x)) or piecewise function
2. **Polygon**: Sympy Polygon object (must be symmetric)
3. **Domain**: Sympy Interval defining the domain of interest (unbounded domains only)

The system then:

1. Finds transition points where the trajectory is parallel to polygon sides
2. Computes active corners for each interval between transition points
3. Determines domain bounds and derivative signs
4. Generates appropriate proof script parameters for unbounded domains

## Functions

### `generate_unbounded_proof_calls(trajectory, poly, domain, x=symbols("x"), y=symbols("y"))`

Generates a list of parameter dictionaries for `unbounded_one_side_proof_script` calls.

**Parameters:**

- `trajectory`: Sympy expression or Piecewise for trajectory (y=f(x))
- `poly`: Sympy Polygon object (must be symmetric)
- `domain`: Sympy Interval domain
- `x, y`: Sympy symbols for variables

**Returns:**

- List of dictionaries containing parameters for proof script calls

**Example:**

```python
from sympy import *
from pvs_utils import generate_unbounded_proof_calls

x, y = symbols("x y")

# Create square polygon
w = 0.5
square_points = [Point(val) for val in [[w, -w], [w, w], [-w, w], [-w, -w]]]
square = Polygon(*square_points)

# Linear trajectory
trajectory = x/2

# Domain (unbounded on left)
domain = Interval(-oo, 5)

# Generate proof calls
proof_calls = generate_unbounded_proof_calls(trajectory, square, domain)

for call in proof_calls:
    print(f"Case label: {call['case_label']}")
    print(f"Deriv lemma: {call['deriv_lemma']}")
    print(f"Max right: {call['max_right']}")
    print(f"Min left: {call['min_left']}")
    print(f"Domain definition: {call['domain_definition']}")
```

### `generate_proof_scripts_from_calls(proof_calls)`

Converts parameter dictionaries into actual proof script strings.

**Parameters:**

- `proof_calls`: List of dictionaries from `generate_unbounded_proof_calls`

**Returns:**

- List of proof script strings

## Generated Parameters

Each proof call dictionary contains:

- **case_label**: String for CASE node (e.g., "max_right >= 5")
- **deriv_lemma**: MVT lemma name (e.g., "mvt_gen_ge_ro_2")
- **max_right**: Maximum right value inside polygon (e.g., "xo + 0.5")
- **min_left**: Minimum left value inside polygon (e.g., "xo - 0.5")
- **domain_definition**: Domain type ("left_open", "right_open", "bounded")
- **max_right_clipped**: Clipped maximum right value
- **min_left_clipped**: Clipped minimum left value
- **interval**: Tuple of interval bounds
- **active_corner**: Active corner for this interval

## Derivative Lemma Naming

The system generates lemma names following the pattern:
`mvt_gen_{deriv_sign}_{domain_type}_2`

Where:

- `deriv_sign`: "ge" (≥) or "le" (≤) based on derivative sign
- `domain_type`: "lo" (left open), "ro" (right open), or "both" (bounded)
- `2`: Always 2 for now

## Domain Types

- **Left open**: `Interval(-oo, b)` unbounded on left
- **Right open**: `Interval(a, oo)` unbounded on right
- **Bounded**: `Interval(a, b)` where both a and b are finite (not currently supported)

## Supported Trajectories

- Functions of x: `y = f(x)`
- Piecewise functions: `Piecewise((expr1, cond1), (expr2, cond2), ...)`
- Both continuous and discontinuous trajectories

## Limitations

- Polygons must be symmetric (central symmetry)
- Trajectories must be functions of x (y=f(x))
- Currently supports only 2D polygons and trajectories
- Derivative evaluation may fail for complex piecewise functions

## Example Usage

See `test_auto_proof_generation.py` for complete examples including:

- Linear trajectories
- Piecewise trajectories
- Unbounded domains
- Different polygon types
