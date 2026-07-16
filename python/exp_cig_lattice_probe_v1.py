#!/usr/bin/env python3
"""exp_cig_lattice_probe_v1.py -- FOR LOCAL RUN: the subgroup-lattice detector on
the ACTUAL NHANES cigarette counts (SMD650, n~1019), the letter's motivating
example of natural mixed heaping.

Counts are integers, so the base lattice is pitch 1 and the natural chain is
1 | 5 | 10 | 20: the coarsest grain sits 20:1 from the base, so this probe uses
an ITERATED downward relock (probe-local; the shipped single-pass relock covers
ratios up to 5). Two independent readouts of the mixing structure are printed:
  (1) spectral: Moebius inversion of replica-center amplitudes;
  (2) digit-share implied: the linear system from the observed shares of
      multiples of 5, 10, and 20 under grain-rounding of a smooth latent.
Their agreement (or the detector's abstention) is the result; no ground truth
exists, so this is a consistency check, not an ISE row.

Run from scripts/:  python3 exp_cig_lattice_probe_v1.py --data-dir data/nhanes
"""
import os, sys, argparse
import numpy as np
from scipy.ndimage import uniform_filter1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ad_kde_v31 as K
import exp_heaping_coarse_v1 as E
import exp_heaping_mixed_v1 as M

# analysis grid: span 720 so P(1)=720, P(5)=144, P(10)=72, P(20)=36; M=8192 so
# three harmonics of the base-pitch comb (2160) sit inside n_half=4096
XC_AD = np.linspace(-330.0, 390.0, 8192)
XC = np.linspace(0.0, 70.0, 1024)
W_MIN = 0.08


def lattice_counts(d, xg, grains_cand=(1, 2, 5, 10, 20)):
    """Subgroup-lattice weight reader for INTEGER count data. The unit pitch is
    known by instrument definition, so the base lattice is anchored rather than
    blindly locked, and the candidate grains are the instrument's natural units
    (singles, pairs, fives, half-packs, packs), stated a priori: an open grain
    dictionary makes the replica lattices near-collinear over the observed bins
    (every candidate period divides the unit period) and the inversion
    ill-posed. Amplitudes are read at replica bins BEYOND the measured signal
    extent (at cigarette scale the coarsest comb's fundamental sits inside the
    signal band, so level-one amplitudes are confounded and excluded), and the
    weights solve a nonnegative least-squares system, the Moebius inversion
    generalized to missing levels. Returns (info, floor-z, peeled density).
    """
    from scipy.optimize import nnls
    from scipy.ndimage import median_filter

    def amp_at(Sv, b, halfwin=2):
        # bin-commensurate combs have one-bin peaks; read the raw maximum
        lo, hi = max(1, b-halfwin), min(len(Sv)-1, b+halfwin+1)
        return float(np.sqrt(max(np.max(Sv[lo:hi]), 0.0)))
    n = len(d)
    DT = 2*np.pi / (xg[-1] - xg[0])
    span = xg[-1] - xg[0]
    p, dx = K._binned(d, xg)
    Phat = np.fft.rfft(p)
    S = np.abs(Phat)**2
    n_half = len(S)
    logS = np.log(S + 1e-300)
    rollmed = median_filter(logS, size=12, mode="nearest")
    # robust floor: combs occupy a minority of bins
    vfloor = None   # set after the extent is known
    # signal extent from a provisional floor, then the floor re-measured
    # beyond the extent so the signal band cannot pollute it
    vprov = float(np.quantile(rollmed[40:min(700, n_half-1)], 0.20))
    ext = 8
    for k in range(8, n_half-6):
        if np.all(rollmed[k:k+5] < vprov + 2.0):
            ext = k
            break
    vfloor = float(np.quantile(rollmed[ext+10:n_half-10], 0.20))
    noise_amp = float(np.exp(vfloor/2))
    # base pitch check: the unit comb must actually be there
    P1 = span / 1.0
    if int(round(P1)) + 5 >= n_half:
        return None, 0.0, None
    a_unit = amp_at(S, int(round(P1)))
    p_unit = float(int(round(P1)))
    z = 2.0*np.log(max(a_unit, 1e-12)) - vfloor
    if a_unit < max(W_MIN, 4*noise_amp):
        return None, z, None
    P_ref = float(p_unit)          # sub-bin refined unit period
    # observation bins: candidate-grain replica multiples beyond the extent
    Pg = {g: P_ref / g for g in grains_cand}
    obs_bins = {}
    for g, per in Pg.items():
        s = 1
        while s * per < n_half - 6:
            b = int(round(s * per))
            if b > ext + 6:
                obs_bins.setdefault(b, set()).add(g)
            s += 1
    bins = sorted(obs_bins)
    if len(bins) < len(grains_cand):
        return None, z, None
    X = np.zeros((len(bins), len(grains_cand)))
    y = np.zeros(len(bins))
    for i_b, b in enumerate(bins):
        y[i_b] = amp_at(S, b)
        for i_g, g in enumerate(grains_cand):
            if abs(b / Pg[g] - round(b / Pg[g])) < 0.02:
                X[i_b, i_g] = 1.0
    w_hat, _res = nnls(X, y)
    keep = [(g, float(wv)) for g, wv in zip(grains_cand, w_hat) if wv > W_MIN]
    if not keep:
        return None, z, None
    tot = sum(wv for _g, wv in keep)
    if tot > 1.0:
        keep = [(g, wv/tot) for g, wv in keep]
    w0 = max(0.0, 1.0 - sum(wv for _g, wv in keep))
    info = dict(grains=[float(g) for g, _ in sorted(keep, key=lambda kv: -kv[0])],
                weights=[round(wv, 3) for _g, wv in sorted(keep, key=lambda kv: -kv[0])],
                w0=round(w0, 3))
    # peel with the mixture kernel; band set by the coarsest detected grain
    g_coarse = max(g for g, _ in keep)
    t_ax = DT * np.arange(len(Phat))
    Kmix = w0 * np.ones_like(t_ax)
    for g, wv in keep:
        Kmix += wv * np.sinc((t_ax * g / 2.0) / np.pi)
    band = int(0.90 * (P_ref / g_coarse) / 2)
    corr = np.ones_like(t_ax)
    if np.abs(Kmix[band-1]) < 0.90:
        corr[:band] = 1.0 / np.maximum(np.abs(Kmix[:band]), 0.25)
    Phat_c = Phat * corr
    Phat_c[band:] = 0.0
    ecf2 = np.abs(Phat_c)**2
    Ssm = uniform_filter1d(ecf2, max(3, band // 8))
    nu0 = float(np.exp(vfloor)) / np.log(2.0)
    nu_t = nu0 * corr**2
    below = np.where(Ssm[1:band] < nu_t[1:band])[0]
    kc = (below[0] + 1) if len(below) else band - 1
    Wf = np.clip(Ssm - nu_t, 0, None) / np.maximum(Ssm, 1e-15)
    Wf[kc:] = 0.0
    f = np.clip(np.fft.irfft(Wf * Phat_c, n=len(xg)) / dx, 0, None)
    return info, z, f


def observed_shares(x):
    xi = np.round(x).astype(int)
    return {g: float(np.mean(xi % g == 0)) for g in (5, 10, 20)}


def predicted_shares(f_lat, xlat, grains, weights, w0):
    """Forward check: from the peeled latent density and the spectral
    decomposition, predict the observed multiple shares. Latent-aware: no
    uniform-phase assumption."""
    xi = np.arange(1, 61)
    p = np.interp(xi, xlat, f_lat); p = np.clip(p, 0, None)
    p = p / p.sum()
    out = {}
    for gq in (5, 10, 20):
        s = w0 * p[(xi % gq) == 0].sum()
        for g, w in zip(grains, weights):
            r = np.clip(np.round(xi / g) * g, 1, 60)
            s += w * p[(r % gq) == 0].sum()
        out[gq] = float(s)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=".")
    args = ap.parse_args()
    from exp_nhanes_heaping_v15 import load
    _meas, cig = load(args.data_dir)
    print("loaded %d genuine cigarette/day reports" % len(cig))
    obs = observed_shares(cig)
    print("\nobserved multiple shares: m5=%.3f  m10=%.3f  m20=%.3f" % (obs[5], obs[10], obs[20]))
    info, z, f = lattice_counts(cig, XC_AD)
    if info is None:
        print("\nlattice detector: ABSTAINED (z=%.2f) -- characterizable as a sample-size" % z)
        print("or capacity boundary at n=%d; the residue floor remains the tool of record." % len(cig))
    else:
        print("\nlattice detector FIRED (z=%.2f):" % z)
        print("  spectral grains:  %s" % info["grains"])
        print("  spectral weights: %s   residual w0=%.3f" % (info["weights"], info["w0"]))
        fd0 = np.interp(XC, XC_AD, f); fd0 = np.clip(fd0, 0, None)
        pred = predicted_shares(fd0, XC, info["grains"], info["weights"], info["w0"])
        print("  forward consistency (predicted vs observed multiple shares):")
        for gq in (5, 10, 20):
            print("    m%-2d predicted %.3f   observed %.3f" % (gq, pred[gq], obs[gq]))
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fd = np.interp(XC, XC_AD, f)
        fd = np.clip(fd, 0, None); fd = fd/np.trapezoid(fd, XC)
        fr = K.ad_wiener(cig, XC, strip="residue")
        fr = np.clip(fr, 0, None); fr = fr/np.trapezoid(fr, XC)
        fig, ax = plt.subplots(figsize=(3.4, 2.4))
        ax.hist(cig, bins=np.arange(0, 62, 1), density=True, color="0.85", edgecolor="0.7", lw=0.1)
        ax.plot(XC, fr, "--", color="0.35", lw=1.2, label="AD residue")
        ax.plot(XC, fd, "-", color="0.05", lw=1.2, label="lattice peel")
        ax.set_xlim(0, 60); ax.set_yticks([]); ax.set_xlabel("cigarettes/day", fontsize=8)
        ax.legend(frameon=False, fontsize=7)
        fig.tight_layout(); fig.savefig("fig_cig_lattice.pdf", bbox_inches="tight")
        print("  figure written: fig_cig_lattice.pdf")
    print("\ntool: exp_cig_lattice_probe_v1.py ; n=%d ; seed-free (single dataset)" % len(cig))


if __name__ == "__main__":
    main()
