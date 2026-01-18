from sympy.plotting import plot
from sympy import symbols, Piecewise

# Define the variable and the piecewise function
x = symbols("x")
trajectory_expr = Piecewise(((x + 1) ** 2, x <= 3), (8 * (x + 1) - 16, x > 3))
trajectory_expr = Piecewise(
    (16.0 - 2.0 * x, (x < 2.0)),
    (18.0 - 2.0 * x, (x < 4.0)),
    (50.0 - 2.0 * x, (x < 21.0)),
    (60.0 - 2.0 * x, (x >= 26.0)),
)
trajectory_expr = Piecewise(
    (14.0 - 1.0 * x, x < 4.0),
    (10.4705882352941 - 0.117647058823529 * x, x < 21.0),
    (15.0 - 0.333333333333333 * x, x >= 21.0),
)
trajectory_expr = Piecewise(
    (14.0 - 1.0 * x, (x >= 0) & (x < 4)),
    (10.48 - 0.12 * x, (x >= 4) & (x < 21)),
    (14.89 - 0.33 * x, (x >= 21) & (x <= 27)),
)


# Plot the piecewise function
p = plot(
    trajectory_expr,
    (x, -2, 30),
    show=True,
    title="Piecewise Trajectory",
    ylabel="y",
    xlabel="x",
)
