#!/usr/bin/env python3
"""
heaping_spectral_v1.py  --  Spectral theory and estimation for heaped data (AD framework)

Purpose
-------
Develop and numerically validate the characteristic-function (CF) theory of heaping
that underlies the AD residue-floor / superposition results in spectral_kde_v3, and
introduce and benchmark a new *AD de-heaping* estimator that treats rounding as
deconvolution against the known box kernel, cut off by the residue-meets-floor rule.

Four studies, all synthetic with closed-form ground truth so ISE is exact:

  V1  Replica/aliasing identity.  Round X to grid D.  The CF of the rounded variable is
        phi_Y(w) = sum_m phi(w + 2*pi*m/D) * sinc((w + 2*pi*m/D) * D/2),  sinc(u)=sin u/u.
      We confirm this to numerical precision, and confirm that dividing the empirical CF
      by sinc(wD/2) recovers phi(w) inside the grid-Nyquist band |w| < pi/D (de-Sheppard),
      while leaving |w| > pi/D non-identified (aliasing).  This is the identifiability limit.

  V2  Estimator benchmark (the headline).  Integrated squared error vs the true density as the
      heaping grid D coarsens, for:
        naive        : Gaussian KDE (Silverman) on the rounded data
        ad_residue   : AD-Wiener with the data-driven residue floor (existing method)
        ad_deheap    : NEW -- residue-floor Wiener applied after box-deconvolution (1/sinc),
                       banded to |w| < pi/D  (this file's proposed estimator)
        deconv_fix   : box-deconvolution with a fixed Silverman cutoff (no residue floor)
                       -- the classical "rounding as uniform measurement error" baseline
        sem_impute   : stochastic-imputation de-heaping (Gross-Rendtel SEM, single known grid)

  V3  Unknown-grid recovery.  Estimate D from the location of the first spectral replica
      (comb tooth) in |phi_hat|^2:  D_hat = 2*pi / w_peak.

  V4  Partial (stochastic) heaping.  Observed law (1-p) f + p f_heaped.  The comb power is
      proportional to p^2; recover p from the replica amplitude and check the estimate.

Conventions (EXPERIMENT-PROTOCOL.md): self-contained (numpy/scipy), fixed seed 20260627,
writes results/heaping_spectral_v1.json, ends with a SUMMARY block to paste back.
Runs fully on synthetic data in-container; no reference software needed for these
floor-and-gain primitives (the matched group of binned data is known to be cyclic).
"""

import json, os, time
import numpy as np

SEED = 20260627
rng = np.random.default_rng(SEED)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------------------
# Ground-truth densities: Gaussian mixtures (closed-form density AND CF).
#   f(x)   = sum_j w_j N(x; mu_j, sig_j^2)
#   phi(t) = sum_j w_j exp(i mu_j t - sig_j^2 t^2 / 2)
# ----------------------------------------------------------------------------
class GM:
    def __init__(self, w, mu, sig, name):
        self.w = np.asarray(w, float); self.w /= self.w.sum()
        self.mu = np.asarray(mu, float); self.sig = np.asarray(sig, float)
        self.name = name
    def pdf(self, x):
        x = np.asarray(x, float)[..., None]
        z = (x - self.mu) / self.sig
        return (self.w / (np.sqrt(2*np.pi)*self.sig) * np.exp(-0.5*z*z)).sum(-1)
    def cf(self, t):
        t = np.asarray(t, float)[..., None]
        return (self.w * np.exp(1j*self.mu*t - 0.5*(self.sig**2)*t*t)).sum(-1)
    def sample(self, n, rng):
        k = rng.choice(len(self.w), size=n, p=self.w)
        return rng.normal(self.mu[k], self.sig[k])

DENSITIES = {
    "gaussian":  GM([1.0], [0.0], [1.0], "Gaussian"),
    "bimodal":   GM([0.5, 0.5], [-1.2, 1.2], [0.5, 0.5], "Bimodal"),
    "kurtotic":  GM([2/3, 1/3], [0.0, 0.0], [1.0, 0.1], "Kurtotic"),   # sharp spike on broad base
    "skewed":    GM([0.2,0.2,0.2,0.2,0.2],
                    [0.0,0.5,1.0833,1.4167,1.6875],
                    [1.0,2/3,4/9,8/27,16/81], "Strongly skewed"),      # Marron-Wand #8
}

# ----------------------------------------------------------------------------
# Grid / FFT plumbing.  Bin samples onto M points; |phi_hat(t_k)| = |FFT(p)_k|
# (Theorem 1 of the paper: eigenvalues of the cyclic group-averaged covariance
#  are the squared ECF).  All spectral filters are real, even functions of w.
# ----------------------------------------------------------------------------
def make_grid(lo, hi, M=2048):
    xg = np.linspace(lo, hi, M, endpoint=False)
    dx = xg[1]-xg[0]
    w  = 2*np.pi*np.fft.fftfreq(M, d=dx)   # angular frequencies, incl. negative
    return xg, dx, w

def bin_prob(samples, xg, dx):
    M = len(xg)
    idx = np.floor((samples - xg[0]) / dx + 0.5).astype(int)
    idx = np.clip(idx, 0, M-1)
    p = np.bincount(idx, minlength=M).astype(float)
    p /= p.sum()
    return p

def ecf_power(p):
    """|phi_hat(w_k)|^2 on the FFT grid."""
    P = np.fft.fft(p)
    return P, np.abs(P)**2

def sinc(u):
    return np.sinc(u/np.pi)   # np.sinc(x)=sin(pi x)/(pi x); we want sin(u)/u

def residue_floor(S, w, frac=0.5):
    """Data-driven max-entropy floor b_hat = median(|phi|^2 in noise band)/ln2.
       Noise band = upper `frac` of |w| (structure has decayed there)."""
    wabs = np.abs(w)
    thr = np.quantile(wabs, 1-frac)
    band = wabs >= thr
    return np.median(S[band]) / np.log(2.0)

def support_mask(S, w, b, c=2.0, smooth=9):
    """Effective-support cutoff (paper Sec IV): smooth |phi|^2, find the first frequency
       (increasing |w|) at which it drops into c*b, keep only frequencies below it.
       This is what discards the high-frequency rounding comb; without it the comb teeth
       leak through the Wiener gain and throw spurious spikes."""
    M = len(w)
    order = np.argsort(np.abs(w))
    Ss = np.convolve(S[order], np.ones(smooth)/smooth, mode="same")
    wcut = np.abs(w).max()
    below = Ss < c*b
    # walk outward from origin; first sustained dip into the floor sets the cutoff
    for i in range(2, M):
        if below[i] and below[min(i+1,M-1)]:
            wcut = np.abs(w[order[i]]); break
    return np.abs(w) < wcut

def reconstruct(P, H, dx):
    """Density on grid = inverse transform of H(w) * phi_hat, clipped & renormalized.
       H real even; P = fft(p).  Returns f_hat on xg."""
    f = np.fft.ifft(P * H).real / dx
    f = np.clip(f, 0, None)
    return f

def renorm(f, dx):
    Z = np.trapezoid(f, dx=dx)
    return f / Z if Z > 0 else f

# ---- estimators (all take rounded samples, return density on xg) ----
def est_naive(samples, xg, dx, w, n):
    p = bin_prob(samples, xg, dx)
    P, S = ecf_power(p)
    h = 1.06 * np.std(samples) * n**(-1/5)     # Silverman
    H = np.exp(-0.5*(h*w)**2)
    return renorm(reconstruct(P, H, dx), dx)

def est_ad_residue(samples, xg, dx, w, n):
    p = bin_prob(samples, xg, dx)
    P, S = ecf_power(p)
    b = residue_floor(S, w)
    keep = support_mask(S, w, b)               # truncate beyond effective support (drops the comb)
    Sstrip = np.clip(S - b, 0, None)
    g = np.where(keep, Sstrip / (Sstrip + b), 0.0)
    return renorm(reconstruct(P, g, dx), dx)

def est_ad_floor(samples, xg, dx, w, n):
    """Bare residue-floor Wiener with effective-support cutoff (the existing spectral method,
       applied WITHOUT box-deconvolution).  Kept to document its failure mode: on coarsely
       rounded (effectively discrete) data |phi_hat|^2 is quasi-periodic and never descends to
       a floor, so the median floor is inflated by the comb and the support cutoff collapses."""
    return est_ad_residue(samples, xg, dx, w, n)

def est_ad_deheap(samples, xg, dx, w, n, D):
    """NEW AD de-heaping estimator.  Rounding to grid D is convolution of f with the width-D
       box (CF = sinc(wD/2)) followed by sampling; so de-heaping is deconvolution by that known
       box CF, valid inside the grid-Nyquist band |w| < pi/D (Theorem: identifiable band).
       Divide phi_hat by sinc, Wiener-taper against the sampling floor 1/n as amplified by the
       deconvolution to (1/n)/sinc^2 (the residual-meets-floor rule sets the cutoff automatically,
       no bandwidth is chosen by hand)."""
    p = bin_prob(samples, xg, dx)
    P, S = ecf_power(p)
    wc = np.pi / D                              # grid Nyquist = edge of the identifiable band
    sc = sinc(w * D / 2.0)
    band = (np.abs(w) < wc) & (np.abs(sc) > 1e-2)
    sc2 = np.clip(sc**2, 1e-4, None)
    floor = (1.0/n) / sc2                        # sampling floor amplified by deconvolution
    Sdh = np.where(band, S / sc2, 0.0)          # deconvolved power estimate |phi|^2
    Sstrip = np.clip(Sdh - floor, 0, None)
    g = np.where(band, Sstrip / (Sstrip + floor), 0.0)   # Wiener gain vs the amplified floor
    Hdh = np.where(band, g / sc, 0.0)           # net filter on phi_hat: gain / sinc
    return renorm(reconstruct(P, Hdh, dx), dx)

def est_deconv_fixed(samples, xg, dx, w, n, D):
    """Classical rounding-as-uniform-error deconvolution: divide by sinc, fixed Silverman
       Gaussian cutoff, no residue floor (shows the instability the floor cures)."""
    p = bin_prob(samples, xg, dx)
    P, S = ecf_power(p)
    h = 1.06 * np.std(samples) * n**(-1/5)
    sc = sinc(w * D / 2.0)
    wc = np.pi / D
    band = (np.abs(w) < wc) & (np.abs(sc) > 0.1)  # 0.1 guard: standard deconvolution ridge cap
    H = np.zeros_like(S)
    H[band] = np.exp(-0.5*(h*w[band])**2) / sc[band]
    return renorm(reconstruct(P, H, dx), dx)

def est_sem_impute(samples_rounded, xg, dx, w, n, D, iters=14, burn=5):
    """Gross-Rendtel-style stochastic-imputation de-heaping, single known grid D.
       S-step: impute latent X_i ~ f_hat restricted to the rounding cell [W-D/2, W+D/2).
       M-step: Silverman Gaussian KDE on the imputed X.  Average post-burn densities.
       Vectorized inverse-CDF sampling over the fixed-width cell segment."""
    M = len(xg)
    W = samples_rounded
    X = W + rng.uniform(-D/2, D/2, size=W.shape)     # init: uniform jitter in cell
    cells = max(1, int(round(D/dx)))
    loi = np.clip(np.floor((W - D/2 - xg[0])/dx).astype(int), 0, M-1)
    cols = loi[:, None] + np.arange(cells)[None, :]   # (n, cells) grid indices
    colc = np.clip(cols, 0, M-1)
    Hgauss = None
    acc = np.zeros(M); nkeep = 0
    for it in range(iters):
        h = 1.06 * np.std(X) * n**(-1/5)
        Hgauss = np.exp(-0.5*(h*w)**2)
        p = bin_prob(X, xg, dx)
        f = np.clip(np.fft.ifft(np.fft.fft(p)*Hgauss).real/dx, 1e-12, None)
        seg = f[colc]                                 # (n, cells) density in each cell
        seg = np.where(cols < M, seg, 0.0)
        c = np.cumsum(seg, axis=1)
        c /= c[:, -1:]
        u = rng.random(len(W))
        k = (c < u[:, None]).sum(axis=1)              # inverse-CDF index within cell
        k = np.clip(k, 0, cells-1)
        X = xg[np.clip(loi + k, 0, M-1)] + rng.uniform(0, dx, size=len(W))
        if it >= burn:
            p = bin_prob(X, xg, dx)
            h = 1.06*np.std(X)*n**(-1/5); Hg = np.exp(-0.5*(h*w)**2)
            acc += np.clip(np.fft.ifft(np.fft.fft(p)*Hg).real/dx, 0, None); nkeep += 1
    return renorm(acc / max(nkeep,1), dx)

def ise(fhat, ftrue, dx):
    return float(np.trapezoid((fhat-ftrue)**2, dx=dx))

# ----------------------------------------------------------------------------
def round_to(x, D):
    return D * np.round(x / D)

def run():
    t0 = time.time()
    results = {"seed": SEED, "meta": "heaping_spectral_v1"}

    # ---------- V1: replica/aliasing identity + de-Sheppard recovery ----------
    gm = DENSITIES["bimodal"]
    D = 0.7
    xg, dx, w = make_grid(-8, 8, 4096)
    # exact CF of rounded variable via the pmf p_k = P(Y=kD) on a fine grid,
    # compared to the Poisson-sum replica formula.
    ks = np.arange(-40, 41)
    # cell probabilities from the true CDF (numeric): integrate pdf over each cell
    xf = np.linspace(-12, 12, 200000); ff = gm.pdf(xf); dxf = xf[1]-xf[0]
    F = np.concatenate([[0], np.cumsum((ff[1:]+ff[:-1])/2*dxf)])
    def cdf(v): return np.interp(v, xf, F)
    pk = cdf((ks+0.5)*D) - cdf((ks-0.5)*D)
    wtest = np.linspace(0.01, np.pi/D*0.995, 200)   # inside grid-Nyquist band
    phiY = np.array([np.sum(pk*np.exp(1j*wt*ks*D)) for wt in wtest])   # exact CF of Y
    # replica formula (truncate m)
    ms = np.arange(-6,7)
    phiY_formula = np.array([np.sum([gm.cf(wt+2*np.pi*m/D)*sinc((wt+2*np.pi*m/D)*D/2) for m in ms]) for wt in wtest])
    rel_replica = np.max(np.abs(phiY-phiY_formula))/np.max(np.abs(phiY))
    # de-Sheppard: phiY / sinc(wD/2) should approx phi(w) in-band
    phi_true = gm.cf(wtest)
    phi_desh = phiY / sinc(wtest*D/2)
    err_desh = np.abs(phi_desh - phi_true)/np.max(np.abs(phi_true))
    err_unc  = np.abs(phiY - phi_true)/np.max(np.abs(phi_true))
    midband = wtest < 0.6*np.pi/D               # comfortably inside Nyquist: aliasing negligible
    results["V1_replica"] = {
        "D": D,
        "rel_err_replica_formula": rel_replica,               # theory identity check
        "deSheppard_midband_|w|<0.6piD": float(err_desh[midband].max()),
        "uncorrected_midband_|w|<0.6piD": float(err_unc[midband].max()),
        "deSheppard_edge_->piD": float(err_desh[~midband].max()),   # aliasing-limited
    }

    # ---------- V2: estimator benchmark, ISE vs coarsening grid ----------
    n = 4000
    Ds = [0.0, 0.25, 0.5, 1.0, 1.5]     # 0.0 = no heaping (control)
    nseed = 12
    xg, dx, w = make_grid(-10, 10, 2048)
    V2 = {}
    for dname, gm in DENSITIES.items():
        ftrue = gm.pdf(xg)
        tab = {}
        for D in Ds:
            errs = {k: [] for k in ["naive","ad_floor","ad_deheap","deconv_fix","sem_impute"]}
            for s in range(nseed):
                sub = np.random.default_rng(SEED + 1000*s + int(D*100))
                x = gm.sample(n, sub)
                y = x if D == 0.0 else round_to(x, D)
                errs["naive"].append(ise(est_naive(y,xg,dx,w,n), ftrue, dx))
                errs["ad_floor"].append(ise(est_ad_floor(y,xg,dx,w,n), ftrue, dx))
                if D == 0.0:
                    # no heaping: de-heap estimators reduce to the bare floor / naive estimators
                    errs["ad_deheap"].append(errs["ad_floor"][-1])
                    errs["deconv_fix"].append(errs["naive"][-1])
                    errs["sem_impute"].append(errs["ad_floor"][-1])
                else:
                    errs["ad_deheap"].append(ise(est_ad_deheap(y,xg,dx,w,n,D), ftrue, dx))
                    errs["deconv_fix"].append(ise(est_deconv_fixed(y,xg,dx,w,n,D), ftrue, dx))
                    errs["sem_impute"].append(ise(est_sem_impute(y,xg,dx,w,n,D), ftrue, dx))
            tab[str(D)] = {k: float(np.mean(v)) for k,v in errs.items()}
        V2[dname] = tab
    results["V2_benchmark_ISE_x1e3"] = {d:{D:{k:round(v*1e3,4) for k,v in row.items()}
                                           for D,row in tab.items()} for d,tab in V2.items()}

    # ---------- V3: unknown-grid recovery from the comb tooth ----------
    V3 = {}
    gm = DENSITIES["bimodal"]
    xg, dx, w = make_grid(-12, 12, 4096)
    for D in [0.5, 1.0, 2.0]:
        rec = []
        for s in range(8):
            sub = np.random.default_rng(SEED + 7*s + int(D*10))
            y = round_to(gm.sample(n, sub), D)
            p = bin_prob(y, xg, dx); _, S = ecf_power(p)
            wpos = w[:len(w)//2]; Spos = S[:len(w)//2]
            # comb tooth ~ near 2pi/D; search a window around it for the local power max
            target = 2*np.pi/D
            win = (wpos > 0.5*target) & (wpos < 1.6*target)
            if win.sum() < 3: continue
            wpk = wpos[win][np.argmax(Spos[win])]
            rec.append(2*np.pi/wpk)
        V3[str(D)] = {"D_true": D, "D_hat_mean": float(np.mean(rec)),
                      "D_hat_std": float(np.std(rec))}
    results["V3_grid_recovery"] = V3

    # ---------- V4: partial heaping, recover the heaping fraction p ----------
    V4 = {}
    gm = DENSITIES["bimodal"]
    xg, dx, w = make_grid(-12, 12, 4096)
    D = 1.0
    for p_true in [0.2, 0.4, 0.6]:
        rec = []
        for s in range(8):
            sub = np.random.default_rng(SEED + 13*s + int(p_true*100))
            x = gm.sample(n, sub)
            heap = sub.random(n) < p_true
            y = x.copy(); y[heap] = round_to(x[heap], D)
            p = bin_prob(y, xg, dx); _, S = ecf_power(p)
            # comb tooth at w=2pi/D.  Only fraction p is heaped, so the m=-1 replica of the
            # heaped part contributes |p*phi(0)|^2 = p^2 at the tooth, over the floor b:
            #   S(2pi/D) ~ p^2 + b   =>   p_hat = sqrt( S_peak - b ).
            b = residue_floor(S, w)
            target = 2*np.pi/D
            win = (np.abs(w) > 0.85*target) & (np.abs(w) < 1.15*target)
            Speak = S[win].max()
            phat = np.sqrt(max(Speak - b, 0.0))
            rec.append(min(phat,1.0))
        V4[str(p_true)] = {"p_true": p_true, "p_hat_mean": float(np.mean(rec)),
                           "p_hat_std": float(np.std(rec))}
    results["V4_partial_heaping"] = V4

    results["runtime_sec"] = round(time.time()-t0, 1)
    with open(os.path.join(OUT, "heaping_spectral_v1.json"), "w") as fh:
        json.dump(results, fh, indent=2)

    # ---------------- SUMMARY (paste this back) ----------------
    print("="*72)
    print("SUMMARY (paste this back)   tool=heaping_spectral_v1   seed=%d" % SEED)
    print("="*72)
    v1 = results["V1_replica"]
    print("[V1] Replica/aliasing identity  (D=%.2f):" % v1["D"])
    print("     replica formula vs exact CF       rel err = %.2e   (expect ~0)" % v1["rel_err_replica_formula"])
    print("     midband |w|<0.6pi/D:  uncorrected = %.2e  ->  after 1/sinc = %.2e" %
          (v1["uncorrected_midband_|w|<0.6piD"], v1["deSheppard_midband_|w|<0.6piD"]))
    print("     band edge ->pi/D (aliasing-limited, irreducible): %.2e" % v1["deSheppard_edge_->piD"])
    print()
    print("[V2] ISE x1e3 vs true density (n=%d, %d seeds). rows=grid D, cols=estimator." % (n, nseed))
    for dname in DENSITIES:
        print("  %-16s  D:     naive  adFlor  ad_DEH  dcnvFx  semImp" % DENSITIES[dname].name)
        for D in Ds:
            r = V2[dname][str(D)]
            print("       %18s%4.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f" %
                  ("", D, r["naive"]*1e3, r["ad_floor"]*1e3, r["ad_deheap"]*1e3,
                   r["deconv_fix"]*1e3, r["sem_impute"]*1e3))
    print()
    print("[V3] Unknown-grid recovery  D_hat = 2pi/w_peak :")
    for D,v in results["V3_grid_recovery"].items():
        print("     D_true=%s  ->  D_hat = %.3f +/- %.3f" % (D, v["D_hat_mean"], v["D_hat_std"]))
    print()
    print("[V4] Partial-heaping fraction recovery (D=1.0):")
    for p,v in results["V4_partial_heaping"].items():
        print("     p_true=%s  ->  p_hat = %.3f +/- %.3f" % (p, v["p_hat_mean"], v["p_hat_std"]))
    print("="*72)
    print("runtime %.1fs   -> results/heaping_spectral_v1.json" % results["runtime_sec"])

if __name__ == "__main__":
    run()
