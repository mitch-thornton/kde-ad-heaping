#!/usr/bin/env python3
"""plot_heaping_v1.py -- figures for the Scientific Reports heaping paper. Seed 20260627.
Writes figures/fig_identifiability.pdf and figures/fig_ise_grid.pdf."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import heaping_spectral_v1 as H
import heaping_spectral_v2 as V2

SEED = 20260627
FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 150})

def fig_identifiability():
    gm = H.DENSITIES["bimodal"]; D = 0.7
    wt = np.linspace(0.01, 1.6*np.pi/D, 400)
    # exact CF of rounded via replica formula
    ms = np.arange(-6, 7)
    phiY = np.array([np.sum([gm.cf(x+2*np.pi*m/D)*H.sinc((x+2*np.pi*m/D)*D/2) for m in ms]) for x in wt])
    phi_true = gm.cf(wt)
    phi_desh = phiY / H.sinc(wt*D/2)
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7))
    ax[0].plot(wt, np.abs(phi_true), "k-", lw=1.6, label=r"$|\varphi|$ true")
    ax[0].plot(wt, np.abs(phiY), "C3--", lw=1.2, label=r"$|\varphi_Y|$ heaped")
    ax[0].plot(wt, np.abs(phi_desh), "C0-.", lw=1.2, label=r"$|\varphi_Y/\mathrm{sinc}|$ de-Sheppard")
    ax[0].axvline(np.pi/D, color="0.4", ls=":", lw=1.0)
    ax[0].text(np.pi/D*1.02, 0.55, r"$\pi/D$", color="0.3")
    ax[0].set_xlabel("frequency $w$"); ax[0].set_ylabel("modulus"); ax[0].set_title("(a) de-Sheppard and the identifiable band")
    ax[0].legend(fontsize=7, loc="upper left")
    # de-Sheppard error vs band fraction
    err = np.abs(phi_desh - phi_true) / np.max(np.abs(phi_true))
    ax[1].semilogy(wt/(np.pi/D), err, "C0-", lw=1.4)
    ax[1].axvline(1.0, color="0.4", ls=":", lw=1.0); ax[1].text(1.03, 1.2e-6, r"$\pi/D$", color="0.3")
    ax[1].set_xlabel(r"$w /(\pi/D)$"); ax[1].set_ylabel("rel. error after 1/sinc")
    ax[1].set_title("(b) exact in-band, aliasing-limited at the edge")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_identifiability.pdf")); plt.close(fig)

def fig_ise_grid():
    gm = H.DENSITIES["bimodal"]; xg, dx, w = H.make_grid(-10, 10, 2048); ftrue = gm.pdf(xg)
    n = 4000; Ds = [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]; nseed = 8
    curves = {k: [] for k in ["naive","ad_deheap","sem","superpose"]}
    for D in Ds:
        e = {k: [] for k in curves}
        for s in range(nseed):
            sub = np.random.default_rng(SEED + 100*s + int(D*100))
            y = H.round_to(gm.sample(n, sub), D)
            e["naive"].append(H.ise(H.est_naive(y,xg,dx,w,n), ftrue, dx))
            e["ad_deheap"].append(H.ise(H.est_ad_deheap(y,xg,dx,w,n,D), ftrue, dx))
            e["sem"].append(H.ise(H.est_sem_impute(y,xg,dx,w,n,D), ftrue, dx))
            e["superpose"].append(H.ise(V2.est_superpose_deheap(y,xg,dx,w,n,D,jitter_rng=np.random.default_rng(3*s+int(D*99))), ftrue, dx))
        for k in curves: curves[k].append(np.mean(e[k]))
    fig, ax = plt.subplots(figsize=(3.6, 2.8))
    sty = {"naive":("C7:","naive KDE"),"ad_deheap":("C0-","AD de-heap"),
           "sem":("C2--","SEM (Gross-Rendtel)"),"superpose":("C3-","superposition de-heap")}
    for k,(s,lab) in sty.items():
        ax.semilogy(Ds, np.array(curves[k])*1e3, s, lw=1.5, label=lab, marker="o", ms=3)
    ax.set_xlabel("heaping grid $D$"); ax.set_ylabel(r"ISE $\times 10^3$ (bimodal)")
    ax.set_title("Estimator ISE vs coarsening")
    ax.legend(fontsize=5.5, loc="upper left", handlelength=1.6, borderpad=0.3, labelspacing=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_ise_grid.pdf")); plt.close(fig)

def fig_spectral_detector():
    """Higher-order comb detection in the spectral basis: symmetry defect vs candidate period."""
    bgi = os.path.abspath(os.path.join(os.path.dirname(__file__),
          "..", "..", "incoming", "1adad087-bgi_v5_16_1_build", "bgi_v5_16_1_build", "reference_software"))
    if not os.path.isdir(bgi):
        print("  (BGI reference software absent; skipping spectral-detector figure)"); return
    sys.path.insert(0, bgi)
    try:
        from bmg.groups import make_quotient, symmetry_defect
    except Exception:
        print("  (BGI import failed; skipping spectral-detector figure)"); return
    M = 512; xg = np.linspace(-12.8, 12.8, M, endpoint=False); dx = xg[1]-xg[0]; n = 4000
    gm = H.DENSITIES["bimodal"]; rng = np.random.default_rng(SEED)
    cand = [8, 16, 24, 32, 48, 64, 96]
    def defects(D):
        x = gm.sample(n, rng); y = x if D is None else D*np.round(x/D)
        idx = np.clip(np.floor((y-xg[0])/dx+0.5).astype(int), 0, M-1)
        p = np.bincount(idx, minlength=M).astype(float); p /= p.sum()
        s = np.abs(np.fft.fft(p))**2; s = s - s.mean()
        R = np.outer(s, s) + 1e-12*np.eye(M)
        return [symmetry_defect(R, make_quotient(M, period=K)) for K in cand]
    fig, ax = plt.subplots(figsize=(3.6, 2.8))
    for D, lab, sty in [(0.8, r"heaped $D{=}0.8$ (true $K{=}32$)", "C0-o"),
                        (1.6, r"heaped $D{=}1.6$ (true $K{=}16$)", "C3-s"),
                        (None, "exact (no heaping)", "C7--^")]:
        ax.plot(cand, defects(D), sty, lw=1.5, ms=4, label=lab)
    ax.set_xlabel("candidate replica period $K$ (spectral bins)")
    ax.set_ylabel(r"symmetry defect $\epsilon_G$")
    ax.set_title("Comb detection in the spectral basis")
    ax.set_ylim(-0.03, 1.12)
    ax.legend(fontsize=5.5, loc="center right", handlelength=1.6, borderpad=0.3, labelspacing=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_spectral_detector.pdf")); plt.close(fig)

if __name__ == "__main__":
    fig_identifiability(); fig_ise_grid(); fig_spectral_detector()
    print("figures written to", os.path.abspath(FIG))
