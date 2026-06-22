import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Style
# -----------------------------
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14

# -----------------------------
# Example circle track
# -----------------------------
np.random.seed(0)

R_true = 300.0
Cx_true = 80.0
Cy_true = 100.0

phi_hits = np.linspace(210, 310, 18) * np.pi / 180.0

x_hits = Cx_true + R_true * np.cos(phi_hits)
z_hits = Cy_true + R_true * np.sin(phi_hits)

x_hits += np.random.normal(0, 3.0, size=len(x_hits))
z_hits += np.random.normal(0, 3.0, size=len(z_hits))

# -----------------------------
# Fixed momentum slice
# -----------------------------
B = 1.0
c_light = 299.792458

p_fixed = c_light * B * (R_true / 1000.0)
R_fixed = p_fixed / (c_light * B) * 1000.0

# -----------------------------
# Hough scan setting
# -----------------------------
rd_values = np.linspace(-250, 100, 400)

hit_indices = np.linspace(0, len(x_hits) - 1, 7, dtype=int)
colors = plt.cm.tab10(np.linspace(0, 1, len(hit_indices)))

theta_peak = np.rad2deg(np.arctan2(Cy_true, Cx_true))
rd_peak = np.sqrt(Cx_true**2 + Cy_true**2) - R_true

# -----------------------------
# Calculate theta-rd curves
# -----------------------------
solutions = []

for color_id, idx in enumerate(hit_indices):
    x = x_hits[idx]
    z = z_hits[idx]

    theta_list = []
    rd_list = []

    for rd in rd_values:
        A = R_fixed + rd
        if abs(A) < 1e-6:
            continue

        rhs = (x**2 + z**2 + A**2 - R_fixed**2) / (2.0 * A)
        norm = np.sqrt(x**2 + z**2)

        if norm < 1e-6:
            continue

        val = rhs / norm

        if abs(val) > 1.0:
            continue

        alpha = np.arctan2(z, x)

        theta1 = alpha + np.arccos(val)
        theta2 = alpha - np.arccos(val)

        for theta in [theta1, theta2]:
            theta = np.arctan2(np.sin(theta), np.cos(theta))
            theta_deg = np.rad2deg(theta)

            theta_list.append(theta_deg)
            rd_list.append(rd)

    solutions.append((np.array(theta_list), np.array(rd_list), color_id))

# -----------------------------
# Draw figure
# -----------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: xTPC-zTPC plane
circle_phi = np.linspace(0, 2 * np.pi, 600)
circle_x = Cx_true + R_true * np.cos(circle_phi)
circle_z = Cy_true + R_true * np.sin(circle_phi)

axes[0].plot(circle_x, circle_z, color="black", linewidth=1.5)
axes[0].scatter(x_hits, z_hits, s=35, color="gray", alpha=0.35)

for idx, c in zip(hit_indices, colors):
    axes[0].scatter(
        x_hits[idx],
        z_hits[idx],
        s=60,
        color=c,
        edgecolor="black",
        linewidth=0.4
    )

axes[0].set_xlabel(r"$x_\mathrm{TPC}$ [mm]")
axes[0].set_ylabel(r"$y_\mathrm{TPC}$ [mm]")
axes[0].set_title(r"Hits in $x_\mathrm{TPC}$-$y_\mathrm{TPC}$ plane")
axes[0].grid(True, alpha=0.5)
axes[0].set_aspect("equal", adjustable="box")

# Right: theta-rd Hough curves
for theta_list, rd_list, color_id in solutions:
    axes[1].scatter(
        theta_list,
        rd_list,
        s=5,
        alpha=0.35,
        color=colors[color_id]
    )



axes[1].set_xlabel(r"$\theta$ [deg]")
axes[1].set_ylabel(r"$r_d$ [mm]")
axes[1].set_title(r"$\theta$-$r_d$ Hough space at fixed $p$")
axes[1].grid(True, alpha=0.5)
axes[1].legend(frameon=False)

for ax in axes:
    ax.xaxis.set_label_coords(0.92, -0.1)
    ax.yaxis.set_label_coords(-0.12, 0.95)

plt.tight_layout()
plt.savefig("hough_circle_xz_thetard_curve.png", dpi=300, bbox_inches="tight")
#plt.savefig("hough_circle_xz_thetard_curve.pdf", bbox_inches="tight")
plt.show()

print(f"fixed p = {p_fixed:.2f} MeV/c")
print(f"maximum: theta = {theta_peak:.2f} deg, rd = {rd_peak:.2f} mm")
