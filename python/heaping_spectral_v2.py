#!/usr/bin/env python3
"""
heaping_spectral_v2.py -- Superposition de-heaping: close the coarse-multimodal gap.

Builds on heaping_spectral_v1 (theory + benchmark) and on Mitch's superpose() machinery
(GMM broad-component base + AD-Wiener residual, exp_datagen_v30.py). The idea:

  The identifiability theorem (v1, Cor. 1) says heaped data determine phi only for |w|<pi/D.
  Beyond that band nothing in the data can recover f; a smooth base must carry it. The
  band-limited de-heaping estimator (v1, ad_deheap) is empty there, which is why iterative
  SEM -- whose imputation supplies a smooth base across the aliasing gap -- wins at coarse
  multimodal grids. So: put a smooth GMM base UNDER the band-limited de-heaped residual.

  f_super : CF = phi_base(w) + g(w)*(phi_deheap(w) - phi_base(w))   for |w| < pi/D,
                 phi_base(w)                                        for |w| >= pi/D.

  phi_base is a BIC Gaussian mixture fit to jitter-dequantized data (y + U(-D/2,D/2)),
  Sheppard-debiased (component variances reduced by D^2/12, since D is known). Beyond the
  identifiable band only the base carries -> it fills exactly the aliasing gap. g(w) is a
  residual Wiener gain that keeps coherent in-band detail above the sampling floor.

Runs fully in-container; seed 20260627.
"""
import json, os, time
import numpy as np
from sklearn.mixture import GaussianMixture
import heaping_spectral_v1 as H

SEED = 20260627
rng = np.random.default_rng(SEED)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT, exist_ok=True)

def gmm_broad_base(y, D, xg, cwidth=1.5, kmax=6, jitter_rng=None):
    """BIC Gaussian mixture on jitter-dequantized data; keep only BROAD components
       (sigma >= cwidth * Silverman h) as the smooth base, mirroring Mitch's superpose().
       Jitter dequantizes (fills the comb); broad-only keeps the base pitch-free and smooth,
       so it carries the aliasing gap |w|>pi/D. No per-component Sheppard de-bias (unstable)."""
    jr = jitter_rng if jitter_rng is not None else np.random.default_rng(0)
    x = (y + jr.uniform(-D/2, D/2, size=len(y)))
    xr = x.reshape(-1, 1)
    best = None
    for k in range(1, kmax+1):
        gm = GaussianMixture(k, covariance_type="full", random_state=0, reg_covar=1e-6).fit(xr)
        b = gm.bic(xr)
        if best is None or b < best[0]:
            best = (b, gm)
    gm = best[1]
    pi = gm.weights_.ravel(); mu = gm.means_.ravel(); var = gm.covariances_.ravel()
    hs = 0.9*min(np.std(x, ddof=1), (np.subtract(*np.percentile(x,[75,25]))/1.349))*max(len(x),2)**(-1/5)
    broad = np.sqrt(var) >= cwidth*hs
    if not broad.any():
        broad = np.array([np.argmax(var)==j for j in range(len(var))])  # keep the broadest
    z = (xg[:,None]-mu[None,:])/np.sqrt(var)[None,:]
    comp = pi[None,:]/np.sqrt(2*np.pi*var)[None,:]*np.exp(-0.5*z*z)
    base = comp[:, broad].sum(1)
    return np.clip(base, 0, None)

def deheap_residual_mass(rm, n, w_r, D):
    """De-heap a residual MASS signal: deconvolve by sinc within |w|<pi/D, Wiener-taper against
       the deconvolution-amplified sampling floor. Mirrors _wiener_filter_mass but with the box
       deconvolution and the grid-Nyquist band. rm real; w_r are nonneg rfft angular freqs."""
    Prm = np.fft.rfft(rm)
    sc = H.sinc(w_r * D / 2.0)
    band = (w_r < np.pi/D) & (np.abs(sc) > 1e-2)
    sc2 = np.clip(sc**2, 1e-4, None)
    Sdh = np.where(band, np.abs(Prm)**2 / sc2, 0.0)
    floor = (1.0/n) / sc2
    gain = np.where(band, np.clip(Sdh-floor,0,None)/(np.clip(Sdh-floor,0,None)+floor), 0.0)
    Hf = np.where(band, gain/sc, 0.0)
    return np.fft.irfft(Hf * Prm, n=len(rm))

def est_superpose_deheap(y, xg, dx, w, n, D, jitter_rng=None):
    p, _dx = H.bin_prob(y, xg, dx), dx
    base = gmm_broad_base(y, D, xg, jitter_rng=jitter_rng)          # smooth base (carries the gap)
    rm = p - base*dx                                                # residual mass
    w_r = 2*np.pi*np.fft.rfftfreq(len(xg), d=dx)
    sharp = deheap_residual_mass(rm, n, w_r, D) / dx                # de-heaped in-band detail
    return H.renorm(np.clip(base + sharp, 0, None), dx)

def run():
    t0 = time.time()
    n = 4000; nseed = 12
    Ds = [0.25, 0.5, 1.0, 1.5]
    xg, dx, w = H.make_grid(-10, 10, 2048)
    dens = ["gaussian","bimodal","kurtotic","skewed"]
    res = {"seed": SEED, "meta": "heaping_spectral_v2_superpose"}
    tab = {}
    for dn in dens:
        gm = H.DENSITIES[dn]; ftrue = gm.pdf(xg); tab[dn] = {}
        for D in Ds:
            e = {k: [] for k in ["naive","ad_deheap","sem_impute","superpose"]}
            for s in range(nseed):
                sub = np.random.default_rng(SEED + 1000*s + int(D*100))
                y = H.round_to(gm.sample(n, sub), D)
                jr = np.random.default_rng(SEED + 5*s + int(D*100))
                e["naive"].append(H.ise(H.est_naive(y,xg,dx,w,n), ftrue, dx))
                e["ad_deheap"].append(H.ise(H.est_ad_deheap(y,xg,dx,w,n,D), ftrue, dx))
                e["sem_impute"].append(H.ise(H.est_sem_impute(y,xg,dx,w,n,D), ftrue, dx))
                e["superpose"].append(H.ise(est_superpose_deheap(y,xg,dx,w,n,D,jitter_rng=jr), ftrue, dx))
            tab[dn][str(D)] = {k: float(np.mean(v)) for k,v in e.items()}
    res["ISE_x1e3"] = {d:{D:{k:round(v*1e3,3) for k,v in row.items()} for D,row in t.items()} for d,t in tab.items()}
    res["runtime_sec"] = round(time.time()-t0,1)
    with open(os.path.join(OUT,"heaping_spectral_v2.json"),"w") as fh: json.dump(res,fh,indent=2)

    print("="*68)
    print("SUMMARY (paste this back)  tool=heaping_spectral_v2  seed=%d" % SEED)
    print("Superposition de-heaping vs naive / ad_deheap / SEM.  ISE x1e3, n=%d, %d seeds." % (n,nseed))
    print("="*68)
    for dn in dens:
        print("  %-15s D:    naive  ad_DEH  semImp  SUPERP" % H.DENSITIES[dn].name)
        for D in Ds:
            r = tab[dn][str(D)]
            best = min(r["ad_deheap"],r["sem_impute"],r["superpose"])
            star = lambda v: "*" if abs(v-best)<1e-12 else " "
            print("       %15s%4.2f %7.2f%s%6.2f%s%6.2f%s%6.2f%s" % ("",D,
                  r["naive"]*1e3," ", r["ad_deheap"]*1e3,star(r["ad_deheap"]),
                  r["sem_impute"]*1e3,star(r["sem_impute"]), r["superpose"]*1e3,star(r["superpose"])))
    print("="*68)
    print("runtime %.1fs  (* = best de-heaper in row)  -> results/heaping_spectral_v2.json" % res["runtime_sec"])

if __name__ == "__main__":
    run()
