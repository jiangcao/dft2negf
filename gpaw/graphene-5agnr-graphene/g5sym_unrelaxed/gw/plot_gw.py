"""g5sym_unrelaxed GW: SCBA convergence and ribbon / lead LDOS, iteration 0 vs last."""
import glob
import json
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ef = json.load(open("../inputs/summary.json"))["efermi"]
dm = np.load("../inputs/device_map.npz")
tag, x = dm["tag"], dm["positions"][:, 0]
e = np.load("outputs/electron_energies.npy") - ef
its = sorted(int(re.findall(r"_(\d+)\.npy", f)[0]) for f in glob.glob("outputs/electron_ldos_*.npy"))
upd = [float(v) for v in re.findall(r"Maximum Self-Energy Update: ([0-9.eE+-]+)", open("gw.out").read())]
ld = {i: np.load(f"outputs/electron_ldos_{i}.npy").mean(1) for i in (its[0], its[-1])}   # ky average
rib = tag == 1
ov = (tag != 1) & (x > dm["positions"][rib, 0].min() - 1) & (x < dm["positions"][rib, 0].max() + 1)

fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw={"width_ratios": [1, 1.6, 1.6], "wspace": 0.28})
ax = axes[0]
ax.semilogy(range(1, len(upd) + 1), upd, marker="o", ms=4, color="#2a78d6", lw=1.2)
ax.axhline(1e-5, color="#52514e", lw=0.7, ls="--")
ax.text(1, 1.4e-5, "convergence_tol 1e-5", fontsize=8, color="#52514e")
ax.set_xlabel("SCBA iteration")
ax.set_ylabel("max self-energy update")
ax.set_title("convergence", fontsize=10)
for a_, sel, title in ((axes[1], rib, "ribbon LDOS (100 LOs)"),
                       (axes[2], ov, "lead LDOS under the ribbon")):
    a_.plot(e, ld[its[0]][:, sel].sum(1), color="#a1a1aa", lw=1.2, label=f"iteration {its[0]} (no GW)")
    a_.plot(e, ld[its[-1]][:, sel].sum(1), color="#eb6834", lw=1.4, label=f"iteration {its[-1]} (GW)")
    a_.axvline(0, color="#52514e", lw=0.7, ls="--")
    a_.set_xlabel(r"$E - E_F$ (eV)")
    a_.set_ylabel("LDOS (arb., k_y average)")
    a_.set_title(title, fontsize=10)
    a_.legend(frameon=False, fontsize=8)
for a_ in axes:
    for s in ("top", "right"):
        a_.spines[s].set_visible(False)
cfg = open("quatrex_config.toml").read()
epsr = re.search(r"^epsilon_r = ([0-9.]+)", cfg, re.M).group(1)
mix = re.search(r"^mixing_factor = ([0-9.]+)", cfg, re.M).group(1)
fig.suptitle(f"g5sym_unrelaxed GW (epsilon_r {epsr}, mixing {mix}, {len(upd)} SCBA iterations, "
             f"last update {upd[-1]:.1e})",
             fontsize=10)
fig.savefig("gw_ldos.png", dpi=160, bbox_inches="tight")
print("wrote gw_ldos.png")
r0, r1 = ld[its[0]][:, rib].sum(1), ld[its[-1]][:, rib].sum(1)
from scipy.signal import find_peaks
for lab, r in (("iter 0", r0), ("last", r1)):
    pk, _ = find_peaks(r, prominence=0.05 * r.max())
    print(lab, "ribbon LDOS peaks (E-EF):", np.round(e[pk], 3).tolist())
