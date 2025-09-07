import numpy as np
import pandas as pd

def relative_from_csv(csv_path, time_grid="as_is"):
    """
    Compute airplane-2's trajectory relative to airplane-1 from a CSV with columns:
      t, f1_lat, f1_lon, f1_alt_ft, f2_lat, f2_lon, f2_alt_ft

    Parameters
    ----------
    csv_path : str
        Path to the CSV file.
    time_grid : {"as_is", "interp_to_f1", "union"}
        - "as_is": assume both aircraft share the same timestamps row-wise (fastest).
        - "interp_to_f1": linearly interpolate aircraft 2 to aircraft 1's t.
        - "union": linearly interpolate both to the union of all unique times.

    Returns
    -------
    pd.DataFrame with columns:
        t               : time
        rel_east_m      : east offset of A2 relative to A1 (meters)
        rel_north_m     : north offset of A2 relative to A1 (meters)
        rel_up_m        : vertical offset of A2 relative to A1 (meters; up positive)
        range_m         : Euclidean separation (meters)
    """
    df = pd.read_csv(csv_path)

    # Basic column sanity check
    required = ["t","f1_lat","f1_lon","f1_alt_ft","f2_lat","f2_lon","f2_alt_ft"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")

    # Ensure numeric types
    for c in required:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=required).sort_values("t")

    # Helper: linear interpolate a (t, series) to new_t (with edge hold)
    def interp_to(new_t, t_src, v_src):
        v = np.interp(new_t, t_src, v_src, left=np.nan, right=np.nan)
        # hold ends if out-of-range
        if np.isnan(v[0]):  v[0]  = v_src.iloc[0]
        if np.isnan(v[-1]): v[-1] = v_src.iloc[-1]
        # fill any internal gaps if any
        if np.isnan(v).any():
            idx = np.flatnonzero(~np.isnan(v))
            v = np.interp(np.arange(len(v)), idx, v[idx])
        return v

    t1 = df["t"].to_numpy()

    if time_grid == "as_is":
        # assume both tracks sampled at the same t in each row
        t = t1
        f1_lat = df["f1_lat"].to_numpy()
        f1_lon = df["f1_lon"].to_numpy()
        f1_alt_ft = df["f1_alt_ft"].to_numpy()
        f2_lat = df["f2_lat"].to_numpy()
        f2_lon = df["f2_lon"].to_numpy()
        f2_alt_ft = df["f2_alt_ft"].to_numpy()

    elif time_grid == "interp_to_f1":
        # interpolate aircraft 2 to aircraft 1's timestamps
        t = t1
        f1_lat = df["f1_lat"].to_numpy()
        f1_lon = df["f1_lon"].to_numpy()
        f1_alt_ft = df["f1_alt_ft"].to_numpy()

        # Build a "traj2 view" grouped by t in case of duplicates
        t2 = df["t"].to_numpy()
        f2_lat = interp_to(t, df["t"], df["f2_lat"])
        f2_lon = interp_to(t, df["t"], df["f2_lon"])
        f2_alt_ft = interp_to(t, df["t"], df["f2_alt_ft"])

    elif time_grid == "union":
        # union of times and interpolate both
        t = np.unique(df["t"].to_numpy())

        for prefix in ("f1_", "f2_"):
            # distinct (t, value) pairs to avoid duplicate-time interpolation issues
            g = df[["t", prefix+"lat", prefix+"lon", prefix+"alt_ft"]].groupby("t").mean().reset_index()
            if prefix == "f1_":
                f1_lat = interp_to(t, g["t"], g[prefix+"lat"])
                f1_lon = interp_to(t, g["t"], g[prefix+"lon"])
                f1_alt_ft = interp_to(t, g["t"], g[prefix+"alt_ft"])
            else:
                f2_lat = interp_to(t, g["t"], g[prefix+"lat"])
                f2_lon = interp_to(t, g["t"], g[prefix+"lon"])
                f2_alt_ft = interp_to(t, g["t"], g[prefix+"alt_ft"])
    else:
        raise ValueError("time_grid must be one of {'as_is','interp_to_f1','union'}")

    # --- Convert per-row lat/lon deltas to local ENU in meters ---
    # Small-angle/equirectangular approximation around A1 at each time.
    # 1 deg latitude ≈ 111_132 meters (sufficient for local separations)
    # 1 deg longitude ≈ 111_320 * cos(phi) meters (phi in radians)
    DEG_TO_RAD = np.pi / 180.0
    m_per_deg_lat = 111_132.0

    phi = f1_lat * DEG_TO_RAD
    m_per_deg_lon = 111_320.0 * np.cos(phi)

    dlat = (f2_lat - f1_lat)  # degrees
    dlon = (f2_lon - f1_lon)  # degrees
    east_m  = dlon * m_per_deg_lon
    north_m = dlat * m_per_deg_lat

    # Vertical (convert feet to meters)
    FT_TO_M = 0.3048
    up_m = (f2_alt_ft - f1_alt_ft) * FT_TO_M

    range_m = np.sqrt(east_m**2 + north_m**2 + up_m**2)

    return pd.DataFrame({
        "t": t,
        "rel_east_m": east_m,
        "rel_north_m": north_m,
        "rel_up_m": up_m,
        "range_m": range_m,
    })

rel = relative_from_csv("solution26_CI0.csv", time_grid="as_is")
rel.to_csv("rel.csv")
print(rel.head())
