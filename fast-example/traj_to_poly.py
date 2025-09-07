import numpy as np
import pandas as pd
import sympy as sp

def _piecewise_linear_sympy(t_vals, y_vals, var_symbol="t", closed_right=True, digits=6,keep_segments=True):
    """
    Build SymPy Piecewise for linear interpolation through (t_i, y_i).
    On [t_i, t_{i+1}): y(t) = y_i + slope_i * (t - t_i)
    The last segment is [t_{n-2}, t_{n-1}] and is closed on the right if closed_right=True.
    """
    df = pd.DataFrame({"t": t_vals, "y": y_vals}).dropna()
    df = df.groupby("t", as_index=False).mean().sort_values("t")
    t = df["t"].to_numpy(float)
    y = df["y"].to_numpy(float)
    if len(t) < 2:
        raise ValueError("Need at least two points to build a piecewise-linear function.")

    T = sp.symbols(var_symbol, real=True)
    def S(x):  
        return sp.Float(float(x), digits)

    pieces = []
    for i in range(len(t)-1):
        t0, t1 = S(t[i]), S(t[i+1])
        y0, y1 = S(y[i]), S(y[i+1])
        slope = (y1 - y0) / (t1 - t0)
        expr  = sp.simplify(y0 + slope*(T - t0))
        if i < len(t)-2 or not closed_right:
            cond = sp.And(T >= t0, T < t1)
        else:
            cond = sp.And(T >= t0, T <= t1) 
        # print(cond)
        pieces.append((expr, cond))
        print(expr)
    # print(len(pieces)) 
    pw = sp.Piecewise(*pieces)
    return pw

def piecewise_from_rel_df(rel_df, axis="east", var_symbol="x", closed_right=True, digits=6):
    """
    rel_df columns expected: t, rel_east_m, rel_north_m, rel_up_m
    axis ∈ {"east","north","up"}
    Returns a SymPy Piecewise for the chosen component.
    """
    colmap = {"east": "rel_east_m", "north": "rel_north_m", "up": "rel_up_m"}
    if axis not in colmap:
        raise ValueError("axis must be one of {'east','north','up'}")
    df = (rel_df[["x", colmap[axis]]]
          .dropna()
          .sort_values("x"))
    return _piecewise_linear_sympy(df["x"].to_numpy(),
                                   df[colmap[axis]].to_numpy(),
                                   var_symbol=var_symbol,
                                   closed_right=closed_right,
                                   digits=digits)


rel = pd.read_csv("rel.csv")

# E_pw = piecewise_from_rel_df(rel, axis="east")
# N_pw = piecewise_from_rel_df(rel, axis="north")
U_pw = piecewise_from_rel_df(rel, axis="up")

print(U_pw)


#U_pw
# Piecewise((0.235477*x - 133.081, (x < 60.0) & (x >= 0)), (0.235477*x - 133.081, ((x >= 60.0) | (x >= 120.0) | (x >= 180.0)) & ((x >= 60.0) | (x >= 120.0) | (x < 240.0)) & ((x >= 60.0) | (x >= 180.0) | (x < 180.0)) & ((x >= 60.0) | (x < 180.0) | (x < 240.0)) & ((x >= 180.0) | (x < 120.0) | (x < 180.0)) & ((x < 120.0) | (x < 180.0) | (x < 240.0))), (0.235477*x - 133.081, (x >= 240.0) & (x < 300.0)), (0.235477*x - 133.081, (x >= 300.0) & (x < 360.0)), (0.235477*x - 133.081, (x >= 360.0) & (x < 420.0)), (0.235477*x - 133.081, ((x >= 420.0) | (x >= 480.0) | (x >= 540.0)) & ((x >= 420.0) | (x >= 480.0) | (x < 600.0)) & ((x >= 420.0) | (x >= 540.0) | (x < 540.0)) & ((x >= 420.0) | (x < 540.0) | (x < 600.0)) & ((x >= 540.0) | (x < 480.0) | (x < 540.0)) & ((x < 480.0) | (x < 540.0) | (x < 600.0))), (0.235477*x - 133.081, (x >= 600.0) & (x < 660.0)), (0.235477*x - 133.081, (x >= 660.0) & (x < 720.0)), (0.235477*x - 133.081, (x >= 720.0) & (x < 780.0)), (0.235477*x - 133.081, (x >= 780.0) & (x < 840.0)), (0.235477*x - 133.081, (x >= 840.0) & (x < 900.0)), (0.235477*x - 133.081, (x >= 900.0) & (x < 960.0)), (0.235477*x - 133.081, (x >= 960.0) & (x < 1020.0)), (0.235477*x - 133.081, (x >= 1020.0) & (x < 1080.0)), (0.235477*x - 133.081, (x >= 1080.0) & (x < 1140.0)), (0.235478*x - 133.082, (x >= 1140.0) & (x <= 1200.0)))