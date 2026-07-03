#!/usr/bin/env python3
"""exp_heaping_mixed_v1.py -- hardened subgroup-lattice peel.
v2 over v1:
  (a) sub-bin refinement of the base period from the ell=1..3 replica centers
      (fixes the ~3% integer-bin grain bias);
  (b) general divisor-lattice scan with Moebius inversion: A_s = sum of w_m over
      m | s, so w_m = sum_{d|m} mu(m/d) A_d; handles non-dyadic nesting;
  (c) downward relock: if significant replica centers exist at P/k for small k,
      the base lattice is re-locked to the coarsest significant spacing (the
      detector may lock a fine grain when a weak coarse comb fails its gates).
Battery adds a non-dyadic nesting case (grains 0.15 | 0.45, factor 3).
"""
import sys, os, json
import numpy as np
from scipy.ndimage import uniform_filter1d
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_heaping_coarse_v1 as E
import ad_kde_v31 as K

SEED0 = 20260627
XG = E.XG
W_MIN = 0.08

def mixed_heap(d, grains, weights, rng):
    gs = list(grains) + [None]; ws = list(weights) + [1.0 - sum(weights)]
    idx = rng.choice(len(gs), size=len(d), p=ws)
    out = d.copy()
    for i, g in enumerate(gs):
        if g is not None:
            m = idx == i; out[m] = np.round(d[m]/g)*g
    return out

def center_amp_pos(S, Q, halfwin=4):
    """(amplitude, refined position) of the replica center near bin Q."""
    Ssm = uniform_filter1d(S, 3)
    lo, hi = max(1, Q-halfwin), min(len(S)-2, Q+halfwin+1)
    j = lo + int(np.argmax(Ssm[lo:hi]))
    y0, y1, y2 = Ssm[j-1], Ssm[j], Ssm[j+1]
    den = (y0 - 2*y1 + y2)
    off = 0.5*(y0-y2)/den if abs(den) > 0 else 0.0
    off = float(np.clip(off, -0.5, 0.5))
    peak = y1 - 0.25*(y0-y2)*off
    return float(np.sqrt(max(peak, 0.0))), j + off

def moebius(n):
    res, m, p = 1, n, 2
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            res = -res
        p += 1
    if m > 1: res = -res
    return res

def ad1_lattice_wiener(d, xg, S_max=12):
    n = len(d)
    DT = 2*np.pi / (xg[-1] - xg[0])
    p, dx = K._binned(d, xg)
    Phat = np.fft.rfft(p)
    S = np.abs(Phat)**2
    P_bins, z, vmed, core = E.detect_period(S)
    if P_bins is None:
        return K.ad_wiener(d, xg, strip="residue"), None, z
    # capacity is deferred: with grain mixtures the base-locked comb may be
    # diluted while another lattice level is strong, so band cleanliness is
    # judged against the strongest lattice amplitude (for single combs the two
    # coincide and the guard is unchanged). Checked after the lattice scan.
    n_half = len(S)
    noise_amp = float(np.exp(vmed/2))
    # ---- (c) downward relock: coarsest significant lattice spacing ----
    # any divisor level with visible amplitude must either be confidently
    # modeled (relock) or force abstention: an unattributed coarse component
    # would alias inside the locked band and corrupt the peel.
    P_base = P_bins
    unmodeled_coarse = False
    for k in (5, 4, 3, 2):
        Pk = P_bins / k
        if Pk < 8: continue
        a1, _ = center_amp_pos(S, int(round(Pk)))
        a2, _ = center_amp_pos(S, int(round(2*Pk)))
        if a1 > max(W_MIN, 6*noise_amp) and a2 > max(W_MIN, 6*noise_amp) and a2 >= 0.8*a1:
            P_base = P_bins / k
            unmodeled_coarse = False
            break
        if a1 > W_MIN and a2 > W_MIN and a2 >= 0.6*a1:
            unmodeled_coarse = True
    if unmodeled_coarse:
        return K.ad_wiener(d, xg, strip="residue"), None, -abs(z)
    # ---- (a) sub-bin refinement from the first replica centers ----
    pos = []
    for ell in (1, 2, 3):
        Q = int(round(ell * P_base))
        if Q + 5 < n_half:
            a, pj = center_amp_pos(S, Q)
            if a > 3*noise_amp:
                pos.append(pj / ell)
    P_ref = float(np.median(pos)) if pos else float(P_base)
    g_base = 2*np.pi / (P_ref * DT)
    # ---- (b) divisor-lattice amplitudes + Moebius inversion ----
    A = {}
    for s in range(1, S_max+1):
        Q = int(round(s * P_ref))
        if Q + 5 >= n_half: break
        A[s], _ = center_amp_pos(S, Q)
    w = {}
    for m in sorted(A):
        wm = sum(moebius(m//dd) * A[dd] for dd in range(1, m+1) if m % dd == 0)
        if wm > W_MIN:
            w[m] = wm
    if not w:
        w = {1: min(A.get(1, 1.0), 1.0)}
    tot = sum(w.values())
    if tot > 1.0:
        w = {m: v/tot for m, v in w.items()}
    w0 = max(0.0, 1.0 - sum(w.values()))
    grains = {m: g_base/m for m in w}
    # band-capacity condition against the strongest lattice amplitude
    A_max = max(max(A.values()), float(np.exp(core/2.0)))
    if 2.0*np.log(max(A_max, 1e-12)) - vmed < 5.0:
        return K.ad_wiener(d, xg, strip="residue"), None, -abs(z)
    # ---- mixture-kernel deconvolution in the coarsest fundamental band ----
    t = DT * np.arange(len(Phat))
    Kmix = w0 * np.ones_like(t)
    for m, wm in w.items():
        Kmix += wm * np.sinc((t * grains[m] / 2.0) / np.pi)
    band = int(0.90 * P_ref / 2)
    corr = np.ones_like(t)
    if np.abs(Kmix[band-1]) < 0.90:
        corr[:band] = 1.0 / np.maximum(np.abs(Kmix[:band]), 0.25)
    Phat_c = Phat * corr
    Phat_c[band:] = 0.0
    ecf2 = np.abs(Phat_c)**2
    Ssm = uniform_filter1d(ecf2, max(3, band // 8))
    nu0 = float(np.exp(vmed)) / np.log(2.0)
    nu_t = nu0 * corr**2
    below = np.where(Ssm[1:band] < nu_t[1:band])[0]
    kc = (below[0] + 1) if len(below) else band - 1
    Wf = np.clip(Ssm - nu_t, 0, None) / np.maximum(Ssm, 1e-15)
    Wf[kc:] = 0.0
    f = np.clip(np.fft.irfft(Wf * Phat_c, n=len(xg)) / dx, 0, None)
    sN = np.trapezoid(f, xg)
    info = dict(grains=[round(grains[m], 4) for m in sorted(w)],
                weights=[round(w[m], 3) for m in sorted(w)], w0=round(w0, 3))
    return (f / sN if sN > 0 else f), info, z

def ise(fh, ft): return float(np.trapezoid((fh-ft)**2, XG))

def run():
    out = {}
    cases = [
        ("mix d .3/.3/.4",   (0.125, 0.25, 0.5), (0.3, 0.3, 0.4)),
        ("mix d .2/.3/.5",   (0.125, 0.25, 0.5), (0.2, 0.3, 0.5)),
        ("mix d+ex .15/.25/.35", (0.125, 0.25, 0.5), (0.15, 0.25, 0.35)),
        ("mix nd .3@.15/.5@.45", (0.15, 0.45), (0.3, 0.5)),
        ("rare-coarse .7/.15/.15", (0.125, 0.25, 0.5), (0.7, 0.15, 0.15)),
        ("single 0.5",       (0.5,), (1.0,)),
        ("exact",            (), ()),
    ]
    print("%-26s %-16s %10s %10s %10s %6s" % ("case","target","residue","single","lattice","fire"))
    for cname, grains, ws in cases:
        for tname in ("Gaussian", "Strongly skewed"):
            comps = E.MW[tname]; ftrue = E.mw_pdf(comps, XG)
            vr, vs, vl, fires, infos = [], [], [], 0, []
            for r in range(15):
                rng = np.random.default_rng(SEED0 + 7919*r + 8000)
                d0 = np.clip(E.mw_sample(comps, 8000, rng), -8+1e-9, 8-1e-9)
                d = mixed_heap(d0, grains, ws, rng) if grains else d0
                vr.append(ise(K.ad_wiener(d, XG, strip="residue"), ftrue))
                fs, _, _ = E.ad1_wiener(d, XG); vs.append(ise(fs, ftrue))
                fl, info, _ = ad1_lattice_wiener(d, XG); vl.append(ise(fl, ftrue))
                fires += (info is not None)
                if info is not None and len(infos) < 1: infos.append(info)
            key = "%s|%s" % (cname, tname)
            out[key] = dict(residue=float(np.mean(vr)), single=float(np.mean(vs)),
                            lattice=float(np.mean(vl)), fires=fires, sample_info=infos)
            print("%-26s %-16s %10.4f %10.4f %10.4f %5d/15" %
                  (cname, tname, 1e3*np.mean(vr), 1e3*np.mean(vs), 1e3*np.mean(vl), fires))
            if infos: print("    e.g. %s" % infos[0])
    json.dump(out, open(os.path.join(HERE, "..", "results", "exp_heaping_mixed_v1.json"), "w"), indent=1)
    print("tool: exp_heaping_mixed_v1.py ; ISE x1e3 ; n=8000 reps=15 seed0=%d" % SEED0)

if __name__ == "__main__":
    run()
