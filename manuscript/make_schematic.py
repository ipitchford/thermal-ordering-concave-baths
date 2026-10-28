"""Schematic (not data) of two relaxation curves that cross: Fig. 0 in Sec. 1.1.
Curves are sums of two exponentials chosen only for illustration."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out = Path(__file__).resolve().parent / "figures"
t = np.linspace(0, 12, 400)
hot = 0.9*np.exp(-2*t) + 0.02*np.exp(-0.3*t)
cold = 0.3*np.exp(-2*t) + 0.08*np.exp(-0.3*t)
fig, ax = plt.subplots(figsize=(4.6, 2.6), layout="constrained")
ax.semilogy(t, hot, color="#7a3b69", label="hotter preparation")
ax.semilogy(t, cold, color="#2a7f74", label="cooler preparation")
i = np.argmin(np.abs(hot - cold)[5:]) + 5
ax.plot(t[i], hot[i], "o", color="black", ms=4)
ax.annotate("crossing", (t[i], hot[i]), xytext=(t[i]+1.2, hot[i]*4), arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=8)
ax.set(xlabel="Time", ylabel="Distance from equilibrium", title="Schematic only: a Mpemba crossing")
ax.legend(frameon=False, fontsize=8)
for ext in ("png", "svg"):
    fig.savefig(out / f"schematic_crossing.{ext}", dpi=200)
