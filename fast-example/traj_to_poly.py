import numpy as np
import pandas as pd
import sympy as sp

def _build_merged_piecewise(t_vals, y_vals, var_symbol="x",
                            closed_right=True, atol=1e-12, rtol=1e-9,
                            digits=None):
    """
    Create a SymPy Piecewise for linear interpolation through (t_i, y_i),
    after MERGING consecutive segments that represent the SAME LINE y = m*t + b
    (within atol/rtol).

    Parameters
    ----------
    t_vals, y_vals : array-like
        Time/value samples. Times must be 1D and (after processing) strictly increasing.
    var_symbol : str
        The symbol to use for the independent variable.
    closed_right : bool
        If True, the last interval is closed on the right: [t_{n-2}, t_{n-1}] with <=.
        Others are [t_i, t_{i+1}) with <.
    atol, rtol : float
        Absolute/relative tolerances for comparing slopes & intercepts.
    digits : int | None
        If not None, convert numeric literals to sp.Float with this many significant digits.
        Keep as None to avoid rounding that might mask small differences.

    Returns
    -------
    sympy.Piecewise
    """
    # Clean + sort + dedup times (average duplicates)
    df = pd.DataFrame({"t": t_vals, "y": y_vals}).dropna()
    df = df.groupby("t", as_index=False).mean().sort_values("t")
    t = df["t"].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=float)
    if len(t) < 2:
        raise ValueError("Need at least two points to build a piecewise-linear function.")

    # Build raw segments: [(t0, t1, m, b)]
    segs = []
    for i in range(len(t)-1):
        dt = t[i+1] - t[i]
        if dt == 0:
            continue
        m = (y[i+1] - y[i]) / dt
        b = y[i] - m * t[i]  # so that y = m*t + b
        segs.append([t[i], t[i+1], m, b])

    if not segs:
        raise ValueError("Degenerate data after deduplication.")

    # Merge consecutive segments with same (m, b) within tolerance
    def close(a, b):
        return abs(a - b) <= (atol + rtol * max(abs(a), abs(b)))

    merged = [segs[0]]
    for t0, t1, m, b in segs[1:]:
        t0_prev, t1_prev, m_prev, b_prev = merged[-1]
        if close(m, m_prev) and close(b, b_prev) and close(t0, t1_prev):
            # Extend the previous segment to the new right endpoint
            merged[-1][1] = t1
        else:
            merged.append([t0, t1, m, b])

    # Build SymPy Piecewise
    T = sp.symbols(var_symbol, real=True)

    def S(x):
        if digits is None:
            return sp.Float(x)
        else:
            return sp.Float(x, digits)

    pieces = []
    for k, (t0, t1, m, b) in enumerate(merged):
        expr = S(m) * T + S(b)
        if (k < len(merged) - 1) or not closed_right:
            cond = sp.And(T >= S(t0), T < S(t1))
        else:
            cond = sp.And(T >= S(t0), T <= S(t1))
        pieces.append((expr, cond))

    # We already merged; let SymPy evaluate normally (no messy unions)
    return sp.Piecewise(*pieces)

# Convenience wrapper for your relative-trajectory DataFrame
def piecewise_from_rel_df_merged(rel_df, axis="up", t_col="x", var_symbol="x",
                                 closed_right=True, atol=1e-12, rtol=1e-9,
                                 digits=None):
    """
    rel_df must contain time column t_col and one of:
      axis='east'  -> 'rel_east_m'
      axis='north' -> 'rel_north_m'
      axis='up'    -> 'rel_up_m'
    Returns a SymPy Piecewise with merged colinear segments.
    """
    colmap = {"east": "relative_lat_ft", "north": "relative_lon_ft", "up": "relative_alt_ft"}
    if axis not in colmap:
        raise ValueError("axis must be one of {'east','north','up'}")

    df = (rel_df[[t_col, colmap[axis]]]
          .dropna()
          .sort_values(t_col)
          .rename(columns={t_col: "t", colmap[axis]: "y"}))

    return _build_merged_piecewise(df["t"].to_numpy(),
                                   df["y"].to_numpy(),
                                   var_symbol=var_symbol,
                                   closed_right=closed_right,
                                   atol=atol, rtol=rtol,
                                   digits=digits)


rel = pd.read_csv("relative_trajectories.csv")  # must have columns like: x, rel_east_m, rel_north_m, rel_up_m
E_pw = piecewise_from_rel_df_merged(rel, axis="east", t_col="t", var_symbol="x",
                                    atol=1e-10, rtol=1e-8, digits=None)
U_pw = piecewise_from_rel_df_merged(rel, axis="up", t_col="t", var_symbol="x",
                                    atol=1e-10, rtol=1e-8, digits=None)
N_pw = piecewise_from_rel_df_merged(rel, axis="north", t_col="t", var_symbol="x",
                                    atol=1e-10, rtol=1e-8, digits=None)
print("Latitude:", E_pw, "\n",  "Longitude:", N_pw, "\n",  "Altitude:", U_pw)
