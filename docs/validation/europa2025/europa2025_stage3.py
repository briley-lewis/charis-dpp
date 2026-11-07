"""Europa Feb 2025: per-cycle disk-integrated Q/U vs time (Mueller on/off) per night, linear trends,
bootstrap of the disk-integrated polarisation, disk-median Q sign; CSV + figure in compare/Europa2025."""
import csv, glob, os, sys, numpy as np
from astropy.io import fits
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
B = os.path.expanduser("~/charis-val/datasets/2025"); OUT = os.path.expanduser("~/charis-val/compare/Europa2025"); os.makedirs(OUT, exist_ok=True)
nights = sys.argv[1:] or ["2025-02-09", "2025-02-10", "2025-02-11"]
yy, xx = np.mgrid[0:201, 0:201]
def disk_mask(I, rad=25):
    w = np.clip(np.nan_to_num(I), 0, None); cx, cy = (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum(); return np.hypot(xx - cx, yy - cy) < rad, (cx, cy)
def per_cycle(d):
    rows = []
    for f in sorted(glob.glob(f"{d}/reduc/proc/pdi_c*_q_derot.fits")):
        c = os.path.basename(f)[5:9]; q = fits.getdata(f, 1); u = fits.getdata(f.replace("_q_", "_u_"), 1)
        iq = fits.getdata(f.replace("_q_", "_iq_"), 1); iu = fits.getdata(f.replace("_q_", "_iu_"), 1); I = 0.5 * (iq + iu)
        h0 = fits.getheader(f, 0); mjd = h0.get("MJD", np.nan)
        Ic = np.nanmean(I, axis=0); m, _ = disk_mask(Ic)
        Q = np.nansum(np.nanmean(q, axis=0)[m]); U = np.nansum(np.nanmean(u, axis=0)[m]); It = np.nansum(Ic[m])
        rows.append((int(c), mjd, Q / It * 100, U / It * 100, np.nanmedian(np.nanmean(q, axis=0)[m])))
    return np.array(rows)
rng = np.random.default_rng(1)
allrows = []
fig, ax = plt.subplots(2, len(nights), figsize=(5 * len(nights), 6), sharey="row", squeeze=False)
print("run                        n   q% mean±sd      u% mean±sd      P%     AoLP   | trend per hour: dq/dt  du/dt (sigma)")
for j, night in enumerate(nights):
    for d in ("py", "py_nomm"):
        D = f"{B}/Europa_{night}/{d}"
        if not glob.glob(f"{D}/reduc/proc/pdi_c*_q_derot.fits"): continue
        a = per_cycle(D); label = f"{night} {d}"
        for r in a: allrows.append([night, d] + list(r))
        t = (a[:, 1] - a[0, 1]) * 24; q = a[:, 2]; u = a[:, 3]
        pq = np.polyfit(t, q, 1); pu = np.polyfit(t, u, 1)
        sq = np.std(q - np.polyval(pq, t)) / np.sqrt(np.sum((t - t.mean()) ** 2)); su = np.std(u - np.polyval(pu, t)) / np.sqrt(np.sum((t - t.mean()) ** 2))
        P = np.hypot(q.mean(), u.mean()); aolp = 0.5 * np.degrees(np.arctan2(u.mean(), q.mean()))
        print(f"{label:26s} {len(a):3d}  {q.mean():+.3f}±{q.std():.3f}   {u.mean():+.3f}±{u.std():.3f}   {P:.3f} {aolp:7.2f}  | {pq[0]:+.4f} ({abs(pq[0]/sq):.1f}σ)  {pu[0]:+.4f} ({abs(pu[0]/su):.1f}σ)")
        n = len(a); bs = []
        for _ in range(2000):
            s = a[rng.integers(0, n, n)]; qq = s[:, 2].mean(); uu = s[:, 3].mean(); bs.append((qq, uu, np.hypot(qq, uu), 0.5 * np.degrees(np.arctan2(uu, qq))))
        bs = np.array(bs); print(f"   bootstrap: q {bs[:,0].mean():+.3f}±{bs[:,0].std():.3f}  u {bs[:,1].mean():+.3f}±{bs[:,1].std():.3f}  P {bs[:,2].mean():.3f}±{bs[:,2].std():.3f}  AoLP {bs[:,3].mean():.2f}±{bs[:,3].std():.2f}")
        qc = fits.getdata(f"{D}/reduc/proc/pdi_q_collapsed.fits", 1); uc = fits.getdata(f"{D}/reduc/proc/pdi_u_collapsed.fits", 1); pic = fits.getdata(f"{D}/reduc/proc/pdi_pi_collapsed.fits", 1)
        iq = [fits.getdata(f, 1) for f in sorted(glob.glob(f"{D}/reduc/proc/pdi_c*_iq_derot.fits"))]; I = np.nanmean(np.stack(iq), axis=(0, 1)); m, c = disk_mask(I, 15)
        print(f"   collapsed: I centroid ({c[0]:.2f},{c[1]:.2f}); disk(r<15) median Q {np.nanmedian(qc[m]):+.4f} U {np.nanmedian(uc[m]):+.4f} PI {np.nanmedian(pic[m]):.4f} I {np.nanmedian(I[m]):.1f}")
        ax[0, j].plot(t, q, ".-", label=d); ax[1, j].plot(t, u, ".-", label=d)
    ax[0, j].set_title(night); ax[1, j].set_xlabel("hours since first cycle"); ax[0, j].legend(fontsize=8)
ax[0, 0].set_ylabel("disk-integrated Q/I [%]"); ax[1, 0].set_ylabel("U/I [%]"); fig.tight_layout(); fig.savefig(f"{OUT}/europa2025_qu_vs_time.png", dpi=120)
with open(f"{OUT}/per_cycle_qu_2025.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["night", "run", "cycle", "mjd", "q_pct", "u_pct", "disk_median_q"]); w.writerows(allrows)
print("->", OUT)
