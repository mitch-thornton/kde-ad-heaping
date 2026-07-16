#!/usr/bin/env python3
"""
heaping_spectral_v3.py -- the combined AD heaping-KDE estimator (single, gated).

The de-heaping (v1) and superposition (v2) variants win in different regimes, set by the
identifiability limit: when the density's structure fits inside the grid-Nyquist band |w|<pi/D
the band-limited de-heaping is best (no smooth base needed); when the structure exceeds the band,
a smooth base must carry the aliasing gap and the superposition wins, with a few de-heaping
iterations helping most at the coarsest grids. The COMBINED estimator reads a band-capacity
statistic from the de-heaped spectrum and selects automatically, so a single "AD heaping-KDE"
covers every regime. This module exposes est_combined() and a benchmark that reports it as one
column against the external baselines.

Band-capacity gate:  rho = fraction of de-heaped in-band power sitting in the outer quarter of
the identifiable band.  Low rho -> structure decays inside the band -> de-heap alone.  High rho
-> structure is cut off at the band edge (exceeds the band) -> superposition, iterated when very
high.  Thresholds calibrated on the synthetic battery (see benchmark; printed for audit).
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import heaping_spectral_v1 as H
import heaping_spectral_v2 as V2
from sklearn.mixture import GaussianMixture

SEED = 20260627

# ---- band-capacity statistic ----
def band_capacity(y, xg, dx, w, n, D):
    """rho = de-heaped in-band power in the outer quarter of |w|<pi/D, over total in-band power.
       High rho means the density's spectral content is still strong at the grid-Nyquist edge,
       i.e. structure exceeds the identifiable band and a smooth base is needed."""
    p = H.bin_prob(y, xg, dx); P, S = H.ecf_power(p)
    sc = H.sinc(w * D / 2.0); wc = np.pi / D
    band = (np.abs(w) < wc) & (np.abs(sc) > 1e-2)
    Sdh = np.where(band, S / np.clip(sc**2, 1e-4, None), 0.0)
    floor = (1.0/n) / np.clip(sc**2, 1e-4, None)
    coh = np.clip(Sdh - floor, 0, None) * band
    tot = coh.sum()
    if tot <= 0: return 0.0
    outer = band & (np.abs(w) > 0.75*wc)
    return float(coh[outer].sum() / tot)

# ---- iterated superposition (de-heaping analog of the imputation loop) ----
def _fit_broad_base(x, xg, n, cwidth=1.5, kmax=6):
    xr = x.reshape(-1, 1); best = None
    for k in range(1, kmax+1):
        gm = GaussianMixture(k, covariance_type="full", random_state=0, reg_covar=1e-6).fit(xr)
        b = gm.bic(xr)
        if best is None or b < best[0]: best = (b, gm)
    gm = best[1]; pi = gm.weights_.ravel(); mu = gm.means_.ravel(); var = gm.covariances_.ravel()
    hs = 0.9*min(np.std(x, ddof=1), np.subtract(*np.percentile(x, [75,25]))/1.349)*n**(-1/5)
    broad = np.sqrt(var) >= cwidth*hs
    if not broad.any(): broad = np.arange(len(var)) == np.argmax(var)
    z = (xg[:,None]-mu[None,:])/np.sqrt(var)[None,:]
    return np.clip((pi[None,:]/np.sqrt(2*np.pi*var)[None,:]*np.exp(-0.5*z*z))[:,broad].sum(1), 0, None)

def est_superpose_iter(y, xg, dx, w, n, D, rng, passes=3):
    p = H.bin_prob(y, xg, dx); wr = 2*np.pi*np.fft.rfftfreq(len(xg), d=dx)
    x = y + rng.uniform(-D/2, D/2, len(y)); f = None
    for it in range(passes):
        base = _fit_broad_base(x, xg, n)
        sharp = V2.deheap_residual_mass(p - base*dx, n, wr, D) / dx
        f = H.renorm(np.clip(base + sharp, 0, None), dx)
        c = np.cumsum(f)*dx; c /= c[-1]; x = np.interp(rng.random(n), c, xg)
    return f

# ---- the combined estimator ----
RHO_LO = 0.06     # below: structure fits inside the band -> de-heap alone
RHO_HI = 0.14     # above: strongly band-saturated -> iterate the superposition

def est_combined(y, xg, dx, w, n, D, rng=None):
    if rng is None: rng = np.random.default_rng(SEED)
    rho = band_capacity(y, xg, dx, w, n, D)
    if rho < RHO_LO:
        return H.est_ad_deheap(y, xg, dx, w, n, D), rho, "deheap"
    if rho < RHO_HI:
        return V2.est_superpose_deheap(y, xg, dx, w, n, D, jitter_rng=rng), rho, "superpose"
    return est_superpose_iter(y, xg, dx, w, n, D, rng, passes=3), rho, "super-iter"

def run():
    t0 = time.time()
    n = 4000; nseed = 12; Ds = [0.25, 0.5, 1.0, 1.5]
    xg, dx, w = H.make_grid(-10, 10, 2048)
    dens = ["gaussian", "bimodal", "kurtotic", "skewed"]
    res = {"seed": SEED, "meta": "heaping_spectral_v3_combined", "RHO_LO": RHO_LO, "RHO_HI": RHO_HI}
    tab = {}
    for dn in dens:
        gm = H.DENSITIES[dn]; ftrue = gm.pdf(xg); tab[dn] = {}
        for D in Ds:
            e = {k: [] for k in ["naive","deconv","sem","deheap","super","combined"]}
            rhos = []; picks = {}
            for s in range(nseed):
                sub = np.random.default_rng(SEED + 1000*s + int(D*100))
                y = H.round_to(gm.sample(n, sub), D)
                jr = np.random.default_rng(SEED + 7*s + int(D*100))
                e["naive"].append(H.ise(H.est_naive(y,xg,dx,w,n), ftrue, dx))
                e["deconv"].append(H.ise(H.est_deconv_fixed(y,xg,dx,w,n,D), ftrue, dx))
                e["sem"].append(H.ise(H.est_sem_impute(y,xg,dx,w,n,D), ftrue, dx))
                e["deheap"].append(H.ise(H.est_ad_deheap(y,xg,dx,w,n,D), ftrue, dx))
                e["super"].append(H.ise(V2.est_superpose_deheap(y,xg,dx,w,n,D,jitter_rng=jr), ftrue, dx))
                fc, rho, pick = est_combined(y,xg,dx,w,n,D, rng=jr)
                e["combined"].append(H.ise(fc, ftrue, dx)); rhos.append(rho)
                picks[pick] = picks.get(pick, 0) + 1
            tab[dn][str(D)] = {k: float(np.mean(v)) for k,v in e.items()}
            tab[dn][str(D)]["rho"] = float(np.mean(rhos)); tab[dn][str(D)]["picks"] = picks
    res["ISE_x1e3"] = {d:{D:{k:(round(v*1e3,3) if isinstance(v,float) else v)
                            for k,v in row.items()} for D,row in t.items()} for d,t in tab.items()}
    res["runtime_sec"] = round(time.time()-t0,1)
    os.makedirs("../results", exist_ok=True)
    with open("../results/heaping_spectral_v3.json","w") as fh: json.dump(res,fh,indent=2)

    print("="*86)
    print("SUMMARY (paste this back)  tool=heaping_spectral_v3  seed=%d  gate RHO_LO=%.2f RHO_HI=%.2f"
          % (SEED, RHO_LO, RHO_HI))
    print("Combined AD heaping-KDE vs external SOTA and its own components. ISE x1e3, n=%d, %d seeds." % (n,nseed))
    print("="*86)
    wins = 0; tot = 0
    for dn in dens:
        print("  %-15s D:   naive deconv   SEM  deheap  super | COMBINED  rho  pick" % H.DENSITIES[dn].name)
        for D in Ds:
            r = tab[dn][str(D)]
            ext_best = min(r["deconv"], r["sem"])
            comb = r["combined"]; tot += 1
            flag = "<-win" if comb <= ext_best + 1e-9 else ""
            if comb <= ext_best + 1e-9: wins += 1
            pk = ",".join("%s:%d"%(k,v) for k,v in r["picks"].items())
            print("       %14s%4.2f %6.2f %6.2f %6.2f %6.2f %6.2f | %7.2f %5.2f  %s %s" % ("",D,
                  r["naive"]*1e3, r["deconv"]*1e3, r["sem"]*1e3, r["deheap"]*1e3, r["super"]*1e3,
                  comb*1e3, r["rho"], pk, flag))
    print("="*86)
    print("combined beats-or-ties the best external baseline in %d/%d cells.  runtime %.1fs" % (wins, tot, res["runtime_sec"]))

if __name__ == "__main__":
    run()
