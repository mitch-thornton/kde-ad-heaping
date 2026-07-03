#!/usr/bin/env python3
"""exp_nhanes_heaping_v16.py -- MIXED-GRAIN study on NHANES measured weights
(companion to exp_nhanes_heaping_v15.py, which remains the single-grain tool of
record; this script runs ONLY the mixed study).

A grain mixture is imposed on the measured adult weights: each report is
rounded to 5, 10, or 20 lb with planted weights, or left exact with a planted
unrounded share, mimicking the mixed digit preference of natural self-reports
(cigarette counts mix multiples of five, ten, and twenty). Arms: the residue
floor, the single-pitch peel forced through the mixture, and the
subgroup-lattice peel of exp_heaping_mixed_v1, which also reads back the grain
weights and the unrounded share. ISE against the no-heaping density; the
mixture assignment is re-drawn over reps for stability.

Run from scripts/:  python3 exp_nhanes_heaping_v16.py --data-dir data/nhanes
"""
import os, sys, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ad_kde_v31 as K
import exp_heaping_mixed_v1 as M
import exp_heaping_coarse_v1 as C
from exp_nhanes_heaping_v15 import load, norm, gkde, ise, XW, XW_AD1

SEED0 = 20260627
GRAINS = (5.0, 10.0, 20.0)
# dominant fine grain (most reports round to the nearest 5), matching natural
# digit preference; the coarse lattice is then found by the downward relock
CASES = [
    ("w=.45/.2/.15, exact .2", (0.45, 0.20, 0.15)),
    ("w=.5/.3/.2",             (0.50, 0.30, 0.20)),
]
REPS = 5


def mixed_heap(d, grains, weights, rng):
    gs = list(grains) + [None]
    ws = list(weights) + [1.0 - sum(weights)]
    idx = rng.choice(len(gs), size=len(d), p=ws)
    out = d.copy()
    for i, g in enumerate(gs):
        if g is not None:
            m = idx == i
            out[m] = np.round(d[m] / g) * g
    return out


def lattice_arm(h):
    """Subgroup-lattice peel on the padded grid, grid-native residue fallback."""
    f_pad, info, _z = M.ad1_lattice_wiener(h, XW_AD1)
    if info is None:
        return norm(K.ad_wiener(h, XW, strip="residue"), XW), None
    f = np.interp(XW, XW_AD1, f_pad)
    return norm(f, XW), info


def single_arm(h):
    f_pad, Dhat, _z = C.ad1_wiener(h, XW_AD1)
    if Dhat is None:
        return norm(K.ad_wiener(h, XW, strip="residue"), XW), None
    return norm(np.interp(XW, XW_AD1, f_pad), XW), Dhat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=".")
    args = ap.parse_args()
    meas, _cig = load(args.data_dir)
    print("loaded %d adult measured weights" % len(meas))
    truth = norm(gkde(meas, K.h_silverman(meas), XW), XW)

    print("\nMIXED-GRAIN study -- grains (5, 10, 20) lb imposed on measured weight;")
    print("ISE x1e3 vs no-heaping density; %d re-draws of the mixture assignment.\n" % REPS)
    for cname, ws in CASES:
        vr, vs, vl, fires = [], [], [], 0
        winfo = []
        for r in range(REPS):
            rng = np.random.default_rng(SEED0 + 7919 * r)
            h = mixed_heap(meas, GRAINS, ws, rng)
            vr.append(1e3 * ise(K.ad_wiener(h, XW, strip="residue"), truth, XW))
            fs, _ = single_arm(h)
            vs.append(1e3 * ise(fs, truth, XW))
            fl, info = lattice_arm(h)
            vl.append(1e3 * ise(fl, truth, XW))
            if info is not None:
                fires += 1
                winfo.append(info)
        print("%s:" % cname)
        print("  residue %8.4f   single-peel %8.4f   lattice %8.4f   fire %d/%d" %
              (np.mean(vr), np.mean(vs), np.mean(vl), fires, REPS))
        if winfo:
            gs = winfo[0]["grains"]
            wmat = np.array([w["weights"] for w in winfo])
            w0s = np.array([w["w0"] for w in winfo])
            print("  recovered grains (lb): %s" % [round(g, 2) for g in gs])
            print("  recovered weights:     %s  (mean over fires)" %
                  [round(float(x), 3) for x in wmat.mean(axis=0)])
            print("  recovered exact share: %.3f   (planted %.2f)" %
                  (float(w0s.mean()), 1.0 - sum(ws)))
        planted = {g: w for g, w in zip(GRAINS, ws)}
        print("  planted: %s + exact %.2f\n" % (planted, 1.0 - sum(ws)))
    print("tool: exp_nhanes_heaping_v16.py ; seed0=%d ; reps=%d" % (SEED0, REPS))


if __name__ == "__main__":
    main()
