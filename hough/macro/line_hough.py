import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Style
# -----------------------------
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14

# -----------------------------
# 1. Generate example x-z hits
# -----------------------------
np.random.seed(0)

p0 = 20.0
p2 = 0.35

z_hits = np.linspace(-250, 250, 15)
x_hits = p0 + p2 * z_hits + np.random.normal(0, 5, size=len(z_hits))

# -----------------------------
# 2. Hough transform
# rho = z*cos(theta) + x*sin(theta)
# -----------------------------
theta_deg = np.linspace(0, 180, 361)
theta_rad = np.deg2rad(theta_deg)

rho_curves = []
for z, x in zip(z_hits, x_hits):
    rho = z * np.cos(theta_rad) + x * np.sin(theta_rad)
    rho_curves.append(rho)

# -----------------------------
# 3. Draw figure
# -----------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

colors = plt.cm.tab20(np.linspace(0, 1, len(z_hits)))

# Left: x-z hit distribution
z_line = np.linspace(z_hits.min(), z_hits.max(), 300)
x_line = p0 + p2 * z_line

axes[0].plot(z_line, x_line, color="black", linewidth=1.5, alpha=0.8)

for z, x, c in zip(z_hits, x_hits, colors):
    axes[0].scatter(z, x, color=c, s=45, edgecolor="black", linewidth=0.3)

axes[0].set_title("Hits in x-z plane")
axes[0].set_xlabel("z [mm]")
axes[0].set_ylabel("x [mm]")
axes[0].grid(True, alpha=0.6)

# Right: Hough space
for rho, c in zip(rho_curves, colors):
    axes[1].plot(theta_deg, rho, color=c, linewidth=1.7)

axes[1].set_title(r"Hough space $(\theta,\rho)$")
axes[1].set_xlabel(r"$\theta$ [deg]")
axes[1].set_ylabel(r"$\rho$ [mm]")
axes[1].grid(True, alpha=0.6)

# -----------------------------
# 4. Move axis titles
# -----------------------------
for ax in axes:
    ax.xaxis.set_label_coords(0.95, -0.07)   # x-axis title to right end
    ax.yaxis.set_label_coords(-0.08, 0.95)   # y-axis title to top

plt.tight_layout()
plt.savefig("hough_xz_example.png", dpi=300, bbox_inches="tight")
plt.savefig("hough_xz_example.pdf", bbox_inches="tight")
plt.show()
