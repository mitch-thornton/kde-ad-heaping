#!/usr/bin/env python3
"""exp_heaping_coarse_v1.py -- spot check: does an AD-1-style first-order stage
improve the KDE residual (and hence floor/bandwidth) on heaped data?

The AD-1 pattern from AD-CLEAN (foundation_ad1.py): detect a deterministic,
group-structured component with a studentized gate; peel it before the
second-order stage; abstain to the unchanged pipeline when the gate fails.

KDE instantiation: rounding samples to a grid of pitch Delta makes the empirical
measure supported on Delta*Z, so its ECF is exactly periodic with period
2*pi/Delta -- a deterministic first-order structure in the spectrum. The stage:
  1. DETECT the replica period blind: score candidate periods P by the
     correlation of the signal band S[0:K] with its replica S[P:P+K],
     studentized against a random-lag null; gate at sig_thresh, else abstain.
  2. PEEL the deterministic rounding kernel: divide phi_hat by
     sinc(t*Delta_hat/2) within the fundamental band (Sheppard/deconvolution
     step with the *detected* kernel), cap the correction, cut beyond the band.
  3. Proceed with the shipped residue floor + soft gain + cutoff.
On exact samples the gate must abstain (pipeline unchanged), mirroring the
AD-CLEAN T2/T3 invariants.

Metric: exact ISE against the TRUE (pre-rounding) density.
Seed of record 20260627. Methods compared:
  fixed    : shipped ad_wiener strip="simple"   (known to fail on heaped data)
  residue  : shipped ad_wiener strip="residue"  (current recommended form)
  ad1      : AD-1 stage + residue flow on the corrected spectrum
"""
import json, os, sys
import numpy as np
from scipy.ndimage import uniform_filter1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ad_kde_v31 as K

SEED0 = 20260627

MW = {
    "Gaussian":        [(1.0, 0.0, 1.0)],
    "Strongly skewed": [(1/8, 3*((2/3)**l - 1), (2/3)**l) for l in range(8)],
    "Claw":            [(0.5, 0.0, 1.0)] + [(0.1, l/2 - 1.0, 0.1) for l in range(5)],
}

def mw_pdf(comps, x):
    out = np.zeros_like(x)
    for w, m, s in comps:
        out += w * np.exp(-0.5*((x-m)/s)**2) / (s*np.sqrt(2*np.pi))
    return out

def mw_sample(comps, n, rng):
    ws = np.array([c[0] for c in comps]); ws /= ws.sum()
    idx = rng.choice(len(comps), size=n, p=ws)
    ms = np.array([c[1] for c in comps]); ss = np.array([c[2] for c in comps])
    return rng.normal(ms[idx], ss[idx])

SPAN = (-8.0, 8.0)
M = 8192
XG = np.linspace(SPAN[0], SPAN[1], M)
DT = 2*np.pi/(SPAN[1]-SPAN[0])          # frequency spacing of rfft bins (rad per bin)

TARGETS = {k: MW[k] for k in ("Gaussian", "Strongly skewed", "Claw")}


def detect_period(S, kmin=16, kmax=1200, K_sig=None,
                  corr_thresh=0.90, z_thresh=5.0, z_harm=3.5, n_harm=3):
    """Blind replica-period detection, gate v3 (resurgence above the quietest
    valley + persistence across harmonics).

    For each candidate P: correlate the signal band's log-spectral shape with
    the segment at P; measure the replica CORE level (median of the top bins);
    studentize its resurgence against the QUIETEST nearby valley, chosen as the
    lower-median of (a) the window just before P and (b) the inter-replica gap
    between P and 2P -- the quietest-window rule keeps a long signal tail from
    inflating the null scale, and (b) makes small periods detectable. Gate at
    z_thresh; require the resurgence to PERSIST at non-decaying level across
    n_harm harmonics (rounding replicas are exactly periodic; a density's own
    beats decay with its component-width envelope). Returns (P, z, vmed) or
    (None, z_best, None) on abstention."""
    n_half = len(S)
    K_sig = K_sig or max(12, n_half // 200)
    logS = np.log(S + 1e-300)

    def seg_stats(P, K_eff, basem, nb):
        seg = logS[P:P+K_eff]
        corr = float(basem @ (seg - seg.mean()) /
                     (nb * (np.linalg.norm(seg - seg.mean()) + 1e-300)))
        core = float(np.median(np.sort(seg)[-max(3, K_eff//4):]))
        return corr, core

    # rolling 12-bin median of the log-spectrum, computed once: valley lookup per
    # candidate is then a min over the precomputed array (quietest-stretch rule)
    W = 12
    from scipy.ndimage import median_filter
    rollmed = median_filter(logS, size=W, mode="nearest")

    def valley_stats(P, K_eff):
        # quietest 12-bin stretch anywhere between the signal band and the second
        # replica: the sliding-median minimum finds genuine floor even when the
        # window regions carry beat or replica-tail leakage.
        cands = []
        for pad in (W//2, W//4, 0):
            a_lo, a_hi = max(K_eff, P - 6 - 36) + pad, P - 6 - pad
            if a_hi > a_lo:
                cands.append((a_lo, a_hi))
            b_lo, b_hi = P + K_eff + pad, min(2*P - 4, n_half - 1) - pad
            if b_hi > b_lo:
                cands.append((b_lo, b_hi))
            if cands:
                break
        if not cands:
            return None, None
        best_c, best_m = None, None
        for lo, hi in cands:
            k = lo + int(np.argmin(rollmed[lo:hi]))
            if best_m is None or rollmed[k] < best_m:
                best_m, best_c = float(rollmed[k]), k
        v = logS[best_c - W//2 : best_c + W//2]
        vmed = best_m
        vsc = 1.4826 * float(np.median(np.abs(v - vmed))) + 1e-9
        return vmed, min(max(vsc, 0.5), 1.2)   # scale within the sampling floor's log-fluctuation range

    kmax = min(kmax, n_half - K_sig - 1)
    cands = []
    for P in range(kmin, kmax):
        K_eff = min(K_sig, max(8, P // 2))
        base = logS[:K_eff]
        basem = base - base.mean()
        nb = np.linalg.norm(basem) + 1e-300
        vmed, vsc = valley_stats(P, K_eff)
        if vmed is None:
            continue
        corr, core = seg_stats(P, K_eff, basem, nb)
        z = (core - vmed) / vsc
        if corr > corr_thresh and z > z_thresh:
            ok, levels = True, [core]
            for h in range(2, n_harm + 1):
                Q = h * P
                if Q + K_eff >= n_half:
                    ok = False; break
                ch, lh = seg_stats(Q, K_eff, basem, nb)
                levels.append(lh)
                if ch < corr_thresh - 0.05 or (lh - vmed) / vsc < z_harm:
                    ok = False; break
            if ok and levels[0] - min(levels) < 2.0:
                cands.append((P, corr, z, vmed))
    if not cands:
        return None, 0.0, None, None
    cset = {c[0] for c in cands}
    fundamentals = [c for c in cands
                    if any(abs(2*c[0] - Q) <= max(2, c[0]//32) for Q in cset)]
    pool = fundamentals if fundamentals else cands
    P, corr, z, vmed = min(pool, key=lambda c: c[0])
    # unbiased floor level for the peel and the capacity test: the MEDIAN of the
    # rolling median over the inter-replica gap (the quietest-stretch MINIMUM used
    # for detection is deliberately sensitivity-biased and must not set the floor)
    K_eff = min(K_sig, max(8, P // 2))
    g_lo, g_hi = P + K_eff, min(2*P - 4, n_half - 1)
    if g_hi - g_lo >= 6:
        vfloor = float(np.median(rollmed[g_lo:g_hi]))
    else:
        a_lo, a_hi = max(K_eff, P - 6 - 36), P - 6
        vfloor = float(np.median(rollmed[a_lo:a_hi])) if a_hi - a_lo >= 6 else float(vmed)
    base = logS[:K_eff]
    basem = base - base.mean()
    nb = np.linalg.norm(basem) + 1e-300
    _c, core = seg_stats(P, K_eff, basem, nb)
    return int(P), float(z), vfloor, float(core)


def ad1_wiener(d, xg):
    """AD-1 stage + residue-flow AD-Wiener. Abstains to the shipped residue
    estimator when no replica period is gated. Works on any uniform grid xg;
    the frequency spacing is derived from the grid span."""
    n = len(d)
    DT = 2*np.pi / (xg[-1] - xg[0])
    p, dx = K._binned(d, xg)
    Phat = np.fft.rfft(p)
    S = np.abs(Phat)**2
    P_bins, z, vmed, core = detect_period(S)
    if P_bins is None:
        return K.ad_wiener(d, xg, strip="residue"), None, z
    # band-capacity condition, n-free: recoverability is a deterministic property,
    # so it is judged deterministic-to-deterministic. With clean separation the
    # replica core stands far above the inter-replica gap floor; aliased overlap
    # fills the gap and shrinks the separation. Below 5 nats the band cannot
    # cleanly hold the signal and the stage abstains.
    if core - vmed < 5.0:
        return K.ad_wiener(d, xg, strip="residue"), None, -abs(z)
    # peel the deterministic rounding kernel: phi_Y = phi_X * sinc(t*Delta/2)
    Delta = 2*np.pi / (P_bins * DT)
    t = DT * np.arange(len(Phat))
    arg = t * Delta / 2.0
    sinc = np.sinc(arg / np.pi)                     # np.sinc(x) = sin(pi x)/(pi x)
    band = int(0.90 * P_bins / 2)                    # stay inside the fundamental band
    corr = np.ones_like(sinc)
    if np.abs(sinc[band-1]) < 0.90:                              # only when attenuation matters
        corr[:band] = 1.0 / np.maximum(np.abs(sinc[:band]), 0.25)   # capped deconvolution
    Phat_c = Phat * corr
    Phat_c[band:] = 0.0                              # cut at the fundamental band edge
    # residue flow on the corrected spectrum; the deconvolution colors the noise,
    # so the floor is frequency-dependent: nu(t) = nu * corr(t)^2
    ecf2 = np.abs(Phat_c)**2
    Ssm = uniform_filter1d(ecf2, max(3, band // 8))   # window scaled to the band, not the grid
    nu0 = float(np.exp(vmed)) / np.log(2.0)
    nu_t = nu0 * corr**2
    below = np.where(Ssm[1:band] < nu_t[1:band])[0]
    kc = (below[0] + 1) if len(below) else band - 1
    Wf = np.clip(Ssm - nu_t, 0, None) / np.maximum(Ssm, 1e-15)
    Wf[kc:] = 0.0
    f = np.clip(np.fft.irfft(Wf * Phat_c, n=len(xg)) / dx, 0, None)
    s = np.trapezoid(f, xg)
    return (f / s if s > 0 else f), Delta, z


def ise(fh, ft, xg):
    return float(np.trapezoid((fh - ft)**2, xg))


def main():
    reps = 15
    out = {"seed0": SEED0, "reps": reps, "results": {}}
    print("%-16s %-8s %-7s | %10s %10s %10s | detect" %
          ("target", "regime", "n", "fixed", "residue", "AD-1"))
    for tname, comps in TARGETS.items():
        ftrue = mw_pdf(comps, XG)
        for regime, Delta in (("exact", None), ("heap0.10", 0.10),
                              ("heap0.25", 0.25), ("heap0.50", 0.50)):
            for n in (2000, 8000):
                r_fix, r_res, r_ad1, det = [], [], [], []
                for r in range(reps):
                    rng = np.random.default_rng(SEED0 + 7919*r + n)
                    d = mw_sample(comps, n, rng)
                    if Delta:
                        d = np.round(d / Delta) * Delta
                    d = np.clip(d, SPAN[0]+1e-9, SPAN[1]-1e-9)
                    r_fix.append(ise(K.ad_wiener(d, XG, strip="simple"), ftrue, XG))
                    r_res.append(ise(K.ad_wiener(d, XG, strip="residue"), ftrue, XG))
                    fh, Dhat, z = ad1_wiener(d, XG)
                    r_ad1.append(ise(fh, ftrue, XG))
                    det.append(Dhat if Dhat else -1)
                det = np.array(det)
                hit = float(np.mean(det > 0))
                dstr = ("abstain %.0f%%" % (100*(1-hit))) if Delta is None else \
                       ("hit %.0f%% Dhat=%.3f" % (100*hit, np.median(det[det > 0]) if hit else -1))
                key = "%s|%s|%d" % (tname, regime, n)
                out["results"][key] = {"fixed": float(np.mean(r_fix)),
                                       "residue": float(np.mean(r_res)),
                                       "ad1": float(np.mean(r_ad1)), "detect": dstr}
                print("%-16s %-8s %-7d | %10.4f %10.4f %10.4f | %s" %
                      (tname, regime, n, 1e3*np.mean(r_fix), 1e3*np.mean(r_res),
                       1e3*np.mean(r_ad1), dstr))
    json.dump(out, open(os.path.join(HERE, "..", "results", "exp_heaping_coarse_v1.json"), "w"), indent=1)
    print("tool: exp_heaping_coarse_v1.py ; seed0=%d ; reps=%d ; ISE x1e3 vs TRUE density" %
          (SEED0, reps))


if __name__ == "__main__":
    main()
