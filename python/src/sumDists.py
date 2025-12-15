#!/usr/bin/env python3
"""
Compare pairwise distance summaries against "truth".

Input files: tab/space-delimited with columns:
leaf1, leaf2, mean, sd, median, min, max, lo95, hi95, lo50, hi50

Assumptions:
- Files are sorted so corresponding rows align.
- The second file's "mean" column is the true value; other columns are ignored.

Outputs one header line (starting with '#') and one results line (tab-separated):
R2_of_means    RMSE    frac_CI95_covers    frac_CI50_covers    frac_within_1sd
"""

import sys
import argparse
from math import isnan, sqrt

def parse_args():
    ap = argparse.ArgumentParser(description="Compare pairwise distance summaries to truth.")
    ap.add_argument("estimates", help="File with Bayesian estimates (means, SDs, CIs).")
    ap.add_argument("truth", help="File with true means (use 'mean' column).")
    ap.add_argument("--allow-mismatch-labels", action="store_true",
                    help="Do not error if leaf-name pairs differ across files; just align by row.")
    return ap.parse_args()

def open_maybe(path):
    return sys.stdin if path == "-" else open(path, "r", encoding="utf-8")

def split_fields(line):
    # tolerate tabs or multiple spaces
    return line.strip().split()

def is_header(line):
    s = line.lstrip()
    return (not s) or s.startswith("#") or s.lower().startswith("leaf") or "leaf" in s.lower()

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return float("nan")

def compute_stats(est_path, truth_path, allow_mismatch):
    n = 0
    sum_sq_err = 0.0

    # For R^2 (Pearson r^2) using two-pass numerically stable approach
    est_means = []
    tru_means = []

    ci95_cov = 0
    ci50_cov = 0
    within_1sd = 0

    with open_maybe(est_path) as fe, open_maybe(truth_path) as ft:
        # iterate lines in lockstep
        while True:
            le = fe.readline()
            lt = ft.readline()
            if not le and not lt:
                break
            if not le or not lt:
                raise ValueError("Files have different number of data lines.")

            if is_header(le) and is_header(lt):
                # skip headers/blank/comment lines in either file
                continue
            # parse fields
            ce = split_fields(le)
            ct = split_fields(lt)
            if len(ce) < 11 or len(ct) < 11:
                raise ValueError("Expected at least 11 columns per line in both files.")

            e_l1, e_l2 = ce[0], ce[1]
            t_l1, t_l2 = ct[0], ct[1]
            if not allow_mismatch and (e_l1 != t_l1 or e_l2 != t_l2):
                raise ValueError(f"Mismatched leaf names at row {n+1}: "
                                 f"({e_l1},{e_l2}) vs ({t_l1},{t_l2})")

            e_mean = safe_float(ce[2])
            e_sd   = safe_float(ce[3])
            e_lo95 = safe_float(ce[7]); e_hi95 = safe_float(ce[8])
            e_lo50 = safe_float(ce[9]); e_hi50 = safe_float(ce[10])

            t_mean = safe_float(ct[2])  # truth

            if any(isnan(v) for v in (e_mean, t_mean)):
                # skip rows with missing core values
                continue

            # accumulate for R^2 and RMSE
            est_means.append(e_mean)
            tru_means.append(t_mean)
            err = e_mean - t_mean
            sum_sq_err += err * err
            n += 1

            # coverage checks (swap if bounds are reversed)
            if not isnan(e_lo95) and not isnan(e_hi95):
                lo95, hi95 = (e_lo95, e_hi95) if e_lo95 <= e_hi95 else (e_hi95, e_lo95)
                if lo95 <= t_mean <= hi95:
                    ci95_cov += 1
            if not isnan(e_lo50) and not isnan(e_hi50):
                lo50, hi50 = (e_lo50, e_hi50) if e_lo50 <= e_hi50 else (e_hi50, e_lo50)
                if lo50 <= t_mean <= hi50:
                    ci50_cov += 1

            # within one SD
            if not isnan(e_sd):
                if abs(err) <= e_sd:
                    within_1sd += 1

    if n == 0:
        raise ValueError("No valid data rows found.")

    # RMSE
    rmse = sqrt(sum_sq_err / n)

    # Pearson r^2
    # r = cov(x,y)/sqrt(var(x)var(y)) ; compute via single pass over centered data
    mx = sum(est_means) / n
    my = sum(tru_means) / n
    sxx = 0.0; syy = 0.0; sxy = 0.0
    for x, y in zip(est_means, tru_means):
        dx = x - mx
        dy = y - my
        sxx += dx*dx
        syy += dy*dy
        sxy += dx*dy
    if sxx <= 0 or syy <= 0:
        r2 = float("nan")
    else:
        r = sxy / sqrt(sxx * syy)
        r2 = r * r

    frac_95 = ci95_cov / n
    frac_50 = ci50_cov / n
    frac_1sd = within_1sd / n

    return r2, rmse, frac_95, frac_50, frac_1sd

def main():
    args = parse_args()
    try:
        r2, rmse, f95, f50, f1sd = compute_stats(args.estimates, args.truth, args.allow_mismatch_labels)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # Output
    header = "#R2_of_means\tRMSE\tfrac_CI95_covers\tfrac_CI50_covers\tfrac_within_1sd"
    print(header)
    # 6 decimal places; adjust if you prefer more/less
    print(f"{r2:.6f}\t{rmse:.6f}\t{f95:.6f}\t{f50:.6f}\t{f1sd:.6f}")

if __name__ == "__main__":
    main()
    
