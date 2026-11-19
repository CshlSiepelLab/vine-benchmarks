#!/usr/bin/env python3
"""Reproduce the centering statistics quoted in VINE_DIST_BIAS_ANALYSIS.md from the
committed calibration table (../data/calibration_data.csv).

Prints, per taxon count and method:
  - overall RMSE, R^2, mean signed error (bias) of the posterior-mean distance
  - OLS calibration fit  est = a + b*true   (b<1 => shrinkage of large distances)
  - fraction of MSE removable by a global linear recalibration
  - mean signed relative bias  Sigma(est-true)/Sigma(true)  by true-distance quintile
  - top-decile shrinkage ratio  mean(est)/mean(true)
"""
import os, math, argparse
import numpy as np, pandas as pd

def rmse(x): return math.sqrt(np.mean(x**2))

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(here, "..", "data", "calibration_data.csv"))
    a = ap.parse_args()
    df = pd.read_csv(a.data)
    df["err"] = df["est"] - df["true"]
    order = ["vine", "beast", "mrbayes"]
    sizes = sorted(df["ntaxa"].unique())

    print("="*92)
    print("OVERALL + CALIBRATION FIT (est = a + b*true)")
    print("="*92)
    h = f"{'ntaxa':>6s} {'method':>8s} | {'RMSE':>7s} {'R2':>6s} {'bias':>8s} | {'slope b':>8s} {'intcpt':>7s} | {'%MSE from miscal':>16s}"
    print(h); print("-"*len(h))
    for nt in sizes:
        for m in order:
            s = df[(df.ntaxa==nt)&(df.method==m)]
            if s.empty: continue
            x, y = s["true"].values, s["est"].values
            b, aic = np.polyfit(x, y, 1)
            resid = y-(aic+b*x)
            raw = rmse(s["err"].values); corr = rmse(resid)
            gap = 100*(raw**2-corr**2)/raw**2
            r2 = np.corrcoef(x, y)[0,1]**2
            print(f"{nt:>6d} {m:>8s} | {raw:7.3f} {r2:6.3f} {s['err'].mean():+8.3f} | {b:8.3f} {aic:+7.3f} | {gap:15.1f}%")
        print()

    print("="*92)
    print("RELATIVE BIAS  Sigma(est-true)/Sigma(true)  BY TRUE-DISTANCE QUINTILE")
    print("="*92)
    for nt in sizes:
        s0 = df[df.ntaxa==nt].copy()
        qs = np.unique(np.quantile(s0[s0.method=="vine"]["true"], np.linspace(0,1,6)))
        qs[0]-=1e-9; qs[-1]+=1e-9
        s0["bin"] = pd.cut(s0["true"], bins=qs)
        print(f"\n#### {nt} taxa ####")
        hdr = f"{'bin':>16s} | " + " | ".join(f"{m:>9s}" for m in order)
        print(hdr)
        for b in sorted(s0["bin"].dropna().unique()):
            cells=[]
            for m in order:
                s=s0[(s0.method==m)&(s0.bin==b)]
                cells.append(f"{s['err'].sum()/s['true'].sum():+9.3f}" if len(s) else "   -")
            print(f"{str(b):>16s} | " + " | ".join(f"{c:>9s}" for c in cells))

    print("\n"+"="*92)
    print("TOP-DECILE SHRINKAGE  mean(est)/mean(true)  (longest 10% of true distances)")
    print("="*92)
    for nt in sizes:
        s0 = df[df.ntaxa==nt]
        thr = np.quantile(s0[s0.method=="vine"]["true"], 0.9)
        row=[f"{nt:>4d} taxa (true>{thr:.2f}):"]
        for m in order:
            s=s0[(s0.method==m)&(s0["true"]>thr)]
            row.append(f"{m}={s['est'].mean()/s['true'].mean():.3f}")
        print("  " + "  ".join(row))

if __name__ == "__main__":
    main()
