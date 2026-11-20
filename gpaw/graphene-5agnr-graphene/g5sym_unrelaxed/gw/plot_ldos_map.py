"""g5sym_unrelaxed GW: LDOS(x, E) heat maps, ribbon and lead layers, first vs last SCBA iteration."""
import glob
import json
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ef = json.load(open("../inputs/summary.json"))["efermi"]
dm = np.load("../inputs/device_map.npz")
tag, x = dm["tag"], dm["positions"][:, 0]
bs = dm["block_sizes"]
e = np.load("outputs/electron_energies.npy") - ef
its = sorted(int(re.findall(r"_(\d+)\.npy", f)[0]) for f in glob.glob("outputs/electron_ldos_*.npy"))
rib = tag == 1
# bins of half a lead period, aligned with the block boundaries (wide gaps
# between atom columns), so every lead bin holds the same two atom columns
dx = float(dm["lead_period"]) / 2
x0 = json.load(open("../inputs/summary.json"))["x_range"][0]
edges = x0 + dx * np.arange(int(np.ceil((x.max() - x0) / dx)) + 1)
xc = (edges[1:] + edges[:-1]) / 2
ib = np.digitize(x, edges) - 1


def binned(ld, sel):
    """Sum the k_y-averaged LDOS of the selected LOs per x bin -> (nE, nbins)."""
    out = np.zeros((len(e), len(xc)))
    np.add.at(out.T, ib[sel], ld[:, sel].T)
    return out


maps = {}
for it in (its[0], its[-1]):
    ld = np.load(f"outputs/electron_ldos_{it}.npy").mean(1)        # (nE, norb)
    maps[it] = {"ribbon": binned(ld, rib), "lead": binned(ld, ~rib)}
vmax = max(m[k].max() for m in maps.values() for k in m)
norm = LogNorm(vmax * 1e-3, vmax)
rx0, rx1 = x[rib].min(), x[rib].max()
pad = np.r_[0, np.cumsum(bs)]
core = (x[pad[2]] if len(bs) > 4 else x.min(), x[pad[-3] - 1])  # inner edges of the 2 copied blocks

fig, axes = plt.subplots(2, 2, figsize=(14, 7.5), sharex=True, sharey=True,
                         gridspec_kw={"hspace": 0.12, "wspace": 0.05})
for col, it in enumerate((its[0], its[-1])):
    for row, layer in enumerate(("ribbon", "lead")):
        ax = axes[row, col]
        A = np.where(maps[it][layer] > 0, maps[it][layer], np.nan)
        im = ax.pcolormesh(xc, e, A, cmap="magma", norm=norm, shading="nearest",
                           rasterized=True)
        ax.axhline(0, color="#ffffff", lw=0.7, ls="--", alpha=0.8)
        for xv in (rx0, rx1):
            ax.axvline(xv, color="#a1a1aa", lw=0.6, ls=":")
        for xv in core:
            ax.axvline(xv, color="#2a78d6", lw=0.8, ls="-", alpha=0.7)
        ax.set_facecolor("#000004")
        lab = "no GW" if it == its[0] else "GW"
        ax.set_title(f"{layer} layer, iteration {it} ({lab})", fontsize=10)
        if col == 0:
            ax.set_ylabel(r"$E - E_F$ (eV)")
        if row == 1:
            ax.set_xlabel(r"x ($\AA$)")
cb = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.01)
cb.set_label(rf"LDOS per {dx:.2f} $\AA$ x bin ($k_y$ average)")
cb.outline.set_visible(False)
fig.suptitle("g5sym_unrelaxed LDOS(x, E): grey dotted = ribbon ends; outside the blue lines the "
             "lead blocks are copies (contacts)", fontsize=10)
fig.savefig("ldos_map.png", dpi=160, bbox_inches="tight")
print("wrote ldos_map.png")
