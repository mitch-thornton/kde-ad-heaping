#!/usr/bin/env python3
"""
exp_nhanes_deheap_v1.py -- NHANES real-data heaping study WITH external baselines.

Answers the SPL rejection's R2 (no external baselines) on real survey microdata:
controlled coarsening of NHANES 2017-2018 measured adult weights, comparing the AD
estimators against the dedicated heaped/measurement-error density literature
(SEM imputation a la Gross-Rendtel; box deconvolution; Sheppard), plus blind grid
and heaped-fraction recovery on the naturally heaped cigarette counts.

Run from scripts/:  python3 exp_nhanes_deheap_v1.py --data-dir ../data
Seed 20260627.  Uses heaping_spectral_v1/v2 estimators (same dir).
"""
import os, sys, json, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import heaping_spectral_v1 as H
import heaping_spectral_v2 as V2
import heaping_spectral_v3 as V3

SEED = 20260627
KG_TO_LB = 2.2046226218

def load(data_dir):
    demo = pd.read_sas(os.path.join(data_dir, "DEMO_J.xpt"))[["SEQN", "RIDAGEYR"]]
    bmx = pd.read_sas(os.path.join(data_dir, "BMX_J.xpt"))[["SEQN", "BMXWT"]]
    df = demo.merge(bmx, on="SEQN")
    df = df[(df["RIDAGEYR"] >= 20) & df["BMXWT"].notna()]
    wt = (df["BMXWT"].to_numpy() * KG_TO_LB)
    smq = pd.read_sas(os.path.join(data_dir, "SMQ_J.xpt"))
    cig = smq["SMD650"].to_numpy()
    cig = cig[np.isfinite(cig)]
    cig = cig[(cig >= 1) & (cig <= 90)]
    return wt, cig

def gauss_kde_ref(d, xg, dx):
    h = 1.06*np.std(d)*len(d)**(-1/5)
    w = 2*np.pi*np.fft.fftfreq(len(xg), d=dx)
    p = H.bin_prob(d, xg, dx)
    return H.renorm(H.reconstruct(np.fft.fft(p), np.exp(-0.5*(h*w)**2), dx), dx)

def run(data_dir):
    rng = np.random.default_rng(SEED)
    wt, cig = load(data_dir)
    n = len(wt)
    print("loaded %d adult measured weights, %d cigarette/day reports" % (n, len(cig)))

    lo, hi = wt.min()-40, wt.max()+40
    M = 4096
    xg = np.linspace(lo, hi, M, endpoint=False); dx = xg[1]-xg[0]
    w = 2*np.pi*np.fft.fftfreq(M, d=dx)
    fref = gauss_kde_ref(wt, xg, dx)                      # no-heaping reference
    def ise(fh): return float(np.trapezoid((fh-fref)**2, dx=dx))

    Ds = [5, 10, 20, 30, 40]
    rows = {}
    for D in Ds:
        y = D*np.round(wt/D)
        r = {}
        r["naive"]     = ise(H.est_naive(y, xg, dx, w, n))
        r["ad_residue"]= ise(H.est_ad_residue(y, xg, dx, w, n))
        r["ad_deheap"] = ise(H.est_ad_deheap(y, xg, dx, w, n, D))
        r["superpose"] = ise(V2.est_superpose_deheap(y, xg, dx, w, n, D, jitter_rng=np.random.default_rng(SEED+D)))
        r["deconv"]    = ise(H.est_deconv_fixed(y, xg, dx, w, n, D))
        r["sem"]       = ise(H.est_sem_impute(y, xg, dx, w, n, D))
        fc,_,_         = V3.est_combined(y, xg, dx, w, n, D, rng=np.random.default_rng(SEED+D))
        r["combined"]  = ise(fc)
        # blind grid recovery from the comb
        p = H.bin_prob(y, xg, dx); _, S = H.ecf_power(p)
        wpos = w[:M//2]; Spos = S[:M//2]; target = 2*np.pi/D
        win = (wpos > 0.5*target) & (wpos < 1.6*target)
        r["Dhat"] = float(2*np.pi/wpos[win][np.argmax(Spos[win])]) if win.sum() >= 3 else float("nan")
        rows[D] = r

    print("\nPART 1 -- measured weight rounded to grid D; ISE x1e3 vs no-heaping reference (n=%d):" % n)
    print("%3s %8s %8s %8s %8s %8s %8s %9s %7s" % ("D","naive","ad_resid","deheap","super","deconv","SEM","COMBINED","Dhat"))
    for D in Ds:
        r = rows[D]
        print("%3d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %9.3f %7.2f" % (D, r["naive"]*1e3, r["ad_residue"]*1e3,
              r["ad_deheap"]*1e3, r["superpose"]*1e3, r["deconv"]*1e3, r["sem"]*1e3, r["combined"]*1e3, r["Dhat"]))

    # cigarettes: digit preference + blind grid/fraction from the comb
    lastdig = (np.round(cig).astype(int) % 10)
    share0 = float(np.mean(lastdig == 0)); share_m10 = float(np.mean(np.round(cig).astype(int) % 10 == 0))
    share_m20 = float(np.mean(np.round(cig).astype(int) % 20 == 0))
    Mc = 2048; xc = np.linspace(0, 95, Mc, endpoint=False); dxc = xc[1]-xc[0]
    wc = 2*np.pi*np.fft.fftfreq(Mc, d=dxc)
    pc = H.bin_prob(cig, xc, dxc); _, Sc = H.ecf_power(pc)
    wcp = wc[:Mc//2]; Scp = Sc[:Mc//2]
    # look for a comb tooth near 2pi/10 (grid 10) and 2pi/5
    def tooth(gr):
        t = 2*np.pi/gr; win = (wcp > 0.8*t) & (wcp < 1.2*t)
        return float(wcp[win][np.argmax(Scp[win])]) if win.sum()>=3 else float("nan")
    ghat10 = 2*np.pi/tooth(10); ghat5 = 2*np.pi/tooth(5)
    print("\nPART 2 -- cigarettes/day (n=%d): multiples of 10 = %.0f%%, of 20 = %.0f%%, last-digit-0 = %.0f%% (uniform 10%%)."
          % (len(cig), share_m10*100, share_m20*100, share0*100))
    print("  blind grid tooth near 10: ghat = %.2f ; near 5: ghat = %.2f" % (ghat10, ghat5))

    out = {"seed": SEED, "n_weight": n, "n_cig": int(len(cig)),
           "weight_ISE_x1e3": {str(D): {k: round(v*1e3,4) if k!="Dhat" else round(v,3)
                                        for k,v in rows[D].items()} for D in Ds},
           "cig": {"share_mult10": share_m10, "share_mult20": share_m20, "share_lastdigit0": share0,
                   "ghat_near10": ghat10, "ghat_near5": ghat5}}
    os.makedirs("../results", exist_ok=True)
    with open("../results/exp_nhanes_deheap_v2.json", "w") as fh: json.dump(out, fh, indent=2)
    print("\n-> ../results/exp_nhanes_deheap_v2.json")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--data-dir", default="../data")
    run(ap.parse_args().data_dir)
