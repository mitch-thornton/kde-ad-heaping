#!/usr/bin/env python3
"""
exp_berlin_age_v1.py -- second real-data example: Berlin resident-register age distribution.

The Berlin resident-registration matrix (EWR201512E, 447 planning areas, 3.61M residents;
the dataset used by Gross and Rendtel for the Kernelheaping package) reports each resident's
age only within a known interval, a mixed-grain grouping (1-, 2-, 3-, 4-, 5-, and 15-year
bins). Grouping is the limiting case of heaping: every value in a bin is collapsed to the bin.
This script (i) pools the age-band counts across all areas, (ii) forms a native-resolution
reference density by dequantizing (placing each count uniformly in its bin) and smoothing,
(iii) imposes coarser UNIFORM grids D and recovers the density with the estimators, measuring
ISE against the native reference, and (iv) detects the imposed grain blind from the comb.

This is the interval-grouping counterpart of the NHANES rounding study; the reference is the
native-resolution density rather than a continuous measurement, stated plainly.

Run from scripts/:  python3 exp_berlin_age_v1.py --data-dir <dir with EWR201512E_Matrix.csv>
Seed 20260627.
"""
import os, sys, csv, re, json, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import heaping_spectral_v1 as H
import heaping_spectral_v2 as V2
import heaping_spectral_v3 as V3

SEED = 20260627

def load_bands(data_dir):
    path = os.path.join(data_dir, "EWR201512E_Matrix.csv")
    rows = list(csv.reader(open(path, encoding="latin-1"), delimiter=';'))
    hdr = rows[0]
    cols = [(i, c) for i, c in enumerate(hdr) if re.fullmatch(r'E_E\d\d_\d+', c)]
    def edges(name):
        a, b = re.fullmatch(r'E_E(\d\d)_(\d+)', name).groups(); return int(a), int(b)
    bands = [(i, *edges(c)) for i, c in cols]
    def num(s): return float(s.replace('.', '').replace(',', '.')) if s else 0.0
    tot = np.zeros(len(bands))
    for r in rows[1:]:
        if len(r) < len(hdr): continue
        for k, (i, a, b) in enumerate(bands):
            tot[k] += num(r[i])
    return [(a, b, c) for (i, a, b), c in zip(bands, tot)], len(rows) - 1

def dequantize_sample(bands, n, rng):
    """Draw n ages by choosing a band in proportion to its count, then uniform within it."""
    los = np.array([a for a, b, c in bands], float)
    his = np.array([b for a, b, c in bands], float)
    w = np.array([c for a, b, c in bands], float); w /= w.sum()
    k = rng.choice(len(bands), size=n, p=w)
    return los[k] + rng.random(n) * (his[k] - los[k])

def run(data_dir):
    rng = np.random.default_rng(SEED)
    bands, n_area = load_bands(data_dir)
    widths = sorted(set(int(b - a) for a, b, c in bands))
    total = sum(c for a, b, c in bands)
    print("Berlin resident register: %d planning areas, %d residents, mixed age-grain widths %s yr"
          % (n_area, int(total), widths))

    n = 40000
    x = dequantize_sample(bands, n, rng)                 # native-resolution pseudo-sample
    lo, hi = -25.0, 135.0; M = 2048
    xg = np.linspace(lo, hi, M, endpoint=False); dx = xg[1] - xg[0]
    w = 2*np.pi*np.fft.fftfreq(M, d=dx)
    # native reference density: Gaussian KDE of the dequantized native-grain sample
    h = 1.06*np.std(x)*n**(-1/5)
    fref = H.renorm(H.reconstruct(np.fft.fft(H.bin_prob(x, xg, dx)), np.exp(-0.5*(h*w)**2), dx), dx)
    def ise(fh): return float(np.trapezoid((fh-fref)**2, dx=dx))

    Ds = [5, 10, 15, 20]                                  # impose coarser uniform grouping (years)
    rows = {}
    for D in Ds:
        y = D*np.round(x/D)
        r = {}
        r["naive"]    = ise(H.est_naive(y, xg, dx, w, n))
        r["deconv"]   = ise(H.est_deconv_fixed(y, xg, dx, w, n, D))
        r["sem"]      = ise(H.est_sem_impute(y, xg, dx, w, n, D))
        r["deheap"]   = ise(H.est_ad_deheap(y, xg, dx, w, n, D))
        r["super"]    = ise(V2.est_superpose_deheap(y, xg, dx, w, n, D, jitter_rng=np.random.default_rng(SEED+D)))
        fc, _, _      = V3.est_combined(y, xg, dx, w, n, D, rng=np.random.default_rng(SEED+D))
        r["comb"]     = ise(fc)
        p = H.bin_prob(y, xg, dx); _, S = H.ecf_power(p)
        wpos = w[:M//2]; Spos = S[:M//2]; target = 2*np.pi/D
        win = (wpos > 0.5*target) & (wpos < 1.6*target)
        r["Dhat"] = float(2*np.pi/wpos[win][np.argmax(Spos[win])]) if win.sum() >= 3 else float("nan")
        rows[D] = r

    print("\nAge density recovered under imposed uniform grouping; ISE x1e3 vs native-resolution reference:")
    print("%4s %8s %8s %8s %8s %8s %9s %7s" % ("D(yr)","naive","deconv","SEM","deheap","super","COMBINED","Dhat"))
    for D in Ds:
        r = rows[D]
        print("%4d %8.3f %8.3f %8.3f %8.3f %8.3f %9.3f %7.2f" % (D, r["naive"]*1e3, r["deconv"]*1e3,
              r["sem"]*1e3, r["deheap"]*1e3, r["super"]*1e3, r["comb"]*1e3, r["Dhat"]))

    out = {"seed": SEED, "n_areas": n_area, "n_residents": int(total), "widths": widths,
           "ISE_x1e3": {str(D): {k: (round(v*1e3,4) if k!="Dhat" else round(v,3))
                                 for k,v in rows[D].items()} for D in Ds}}
    os.makedirs("../results", exist_ok=True)
    with open("../results/exp_berlin_age_v1.json", "w") as fh: json.dump(out, fh, indent=2)
    print("\n-> ../results/exp_berlin_age_v1.json")

    # figure: age density recovery at a coarse imposed grain
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
        Dp = 20; y = Dp*np.round(x/Dp)
        f_naive = H.est_naive(y, xg, dx, w, n)
        f_comb, _, _ = V3.est_combined(y, xg, dx, w, n, Dp, rng=np.random.default_rng(SEED+Dp))
        fig, ax = plt.subplots(figsize=(3.7, 2.8))
        m = (xg >= 0) & (xg <= 100)
        ax.plot(xg[m], fref[m], "k-", lw=1.6, label="native-grain reference")
        ax.plot(xg[m], f_naive[m], "C7:", lw=1.4, label="naive KDE (grouped to 20 yr)")
        ax.plot(xg[m], f_comb[m], "C3-", lw=1.4, label="combined de-heaping")
        ax.set_xlabel("age (years)"); ax.set_ylabel("density")
        ax.set_title("Berlin register age density, 20-year grouping")
        ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig("../figures/fig_berlin_age.pdf"); plt.close(fig)
        print("-> ../figures/fig_berlin_age.pdf")
    except Exception as e:
        print("  (figure skipped:", e, ")")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="../../incoming/berlindata/berlin-data")
    run(os.path.abspath(ap.parse_args().data_dir))
