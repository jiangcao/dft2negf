"""T(E) of the g5sym_unrelaxed QTBM run, per transverse ky and the ky average."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

cfg_ef = json.load(open("../inputs/summary.json"))["efermi"]
e = np.load("outputs/electron_energies.npy") - cfg_ef
T = np.asarray(np.load("outputs/transmission_lr.npy")).real    # (nky, nE)
nky = T.shape[0]
ky = (np.arange(nky) + 0.5) / nky - 0.5
fig, ax = plt.subplots(figsize=(7.5, 4.2))
cols = ["#2a78d6", "#eb6834", "#1baf7a"]
for i in range(nky):
    ax.plot(e, T[i], lw=1.0, color=cols[i % 3], label=f"$k_y$ = {ky[i]:+.3f}")
ax.plot(e, T.mean(0), lw=2.0, color="#0b0b0b", label="average")
ax.axvline(0, color="#52514e", lw=0.7, ls="--")
ax.set_yscale("log")
ax.set_ylim(1e-5, None)
ax.set_xlabel(r"$E - E_F$ (eV)")
ax.set_ylabel("transmission per transverse cell")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=8)
ax.set_title("g5sym_unrelaxed: ballistic T(E), LO model (k3x3 correction)", fontsize=10)
fig.savefig("transmission.png", dpi=160, bbox_inches="tight")
print("wrote transmission.png")
