"""Aperture polarisation of a point source from charis pdi products: sums inside r < 12 px of the
collapsed Q and U and of the mean I_Q/I_U double sum; P and AoLP, plus P per channel."""
import glob, os, sys, numpy as np
from astropy.io import fits
yy, xx = np.mgrid[0:201, 0:201]
for d in sys.argv[1:]:
    proc = f"{d}/reduc/proc"
    if not os.path.exists(f"{proc}/pdi_q_collapsed.fits"): print(d, "no products"); continue
    q = fits.getdata(f"{proc}/pdi_q_collapsed.fits", 1); u = fits.getdata(f"{proc}/pdi_u_collapsed.fits", 1)
    qc = fits.getdata(f"{proc}/pdi_q.fits", 1); uc = fits.getdata(f"{proc}/pdi_u.fits", 1)
    iq = np.stack([fits.getdata(f, 1) for f in sorted(glob.glob(f"{proc}/pdi_c*_iq_derot.fits"))])
    iu = np.stack([fits.getdata(f, 1) for f in sorted(glob.glob(f"{proc}/pdi_c*_iu_derot.fits"))])
    Icube = 0.5 * (np.nanmean(iq, axis=0) + np.nanmean(iu, axis=0)); I = np.nanmean(Icube, axis=0)
    w = np.clip(np.nan_to_num(I), 0, None); cx, cy = (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()
    m = np.hypot(xx - cx, yy - cy) < 12
    Q, U, It = np.nansum(q[m]), np.nansum(u[m]), np.nansum(I[m])
    P = np.hypot(Q, U) / It * 100; aolp = 0.5 * np.degrees(np.arctan2(U, Q))
    ncyc = iq.shape[0]
    per = []
    for il in (0, 5, 10, 15, 21):
        Ql, Ul, Il = np.nansum(qc[il][m]), np.nansum(uc[il][m]), np.nansum(Icube[il][m]); per.append(np.hypot(Ql, Ul) / Il * 100)
    print(f"{d:10s} cycles {ncyc:3d}  centre ({cx:.1f},{cy:.1f})  P {P:.3f} %  AoLP {aolp:7.2f} deg  q {Q/It*100:+.3f} u {U/It*100:+.3f}  | P per channel 0/5/10/15/21: " + " ".join(f"{p:.2f}" for p in per))
