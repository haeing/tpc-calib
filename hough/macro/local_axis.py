import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 13

# -----------------------------
# Helix parameters
# -----------------------------
R = 160.0
Cx = 0.0
Cy = 0.0

Z0 = 0.0
Dz = 0.08

phi = np.linspace(np.deg2rad(240), np.deg2rad(300), 100)
phi_ref = np.mean(phi)

x = Cx + R*np.cos(phi)
y = Cy + R*np.sin(phi)

# angle difference around phi_ref, kept near zero
dphi = np.arctan2(np.sin(phi - phi_ref), np.cos(phi - phi_ref))
s = R*dphi
z = Z0 + Dz*s

# -----------------------------
# Hits
# -----------------------------
np.random.seed(1)

hit_idx = np.linspace(0, len(phi)-1, 14, dtype=int)

xh = x[hit_idx] + np.random.normal(0, 2.0, len(hit_idx))
yh = y[hit_idx] + np.random.normal(0, 2.0, len(hit_idx))

phih = np.arctan2(yh - Cy, xh - Cx)
dphih = np.arctan2(np.sin(phih - phi_ref), np.cos(phih - phi_ref))
sh = R*dphih
zh = Z0 + Dz*sh + np.random.normal(0, 0.8, len(hit_idx))

colors = plt.cm.tab20(np.linspace(0, 1, len(hit_idx)))

i_rep = len(hit_idx)//2
x0, y0, z0 = xh[i_rep], yh[i_rep], zh[i_rep]
s0 = sh[i_rep]

# -----------------------------
# Figure layout
# -----------------------------
fig = plt.figure(figsize=(12, 8))
gs = gridspec.GridSpec(
    2, 2,
    width_ratios=[1.0, 1.25],
    height_ratios=[1.0, 0.75],
    wspace=0.25,
    hspace=0.32
)

# -----------------------------
# xTPC-yTPC projection
# -----------------------------
ax_xy = fig.add_subplot(gs[0, 0])

circle_phi = np.linspace(0, 2*np.pi, 500)
circle_x = Cx + R*np.cos(circle_phi)
circle_y = Cy + R*np.sin(circle_phi)

ax_xy.plot(circle_x, circle_y, color="black", lw=1.5)
ax_xy.plot(x, y, color="black", lw=2.0, alpha=0.5)

for xx, yy, c in zip(xh, yh, colors):
    ax_xy.scatter(xx, yy, s=55, color=c, edgecolor="black", linewidth=0.5)

ax_xy.scatter(Cx, Cy, marker="x", s=100, color="black")
ax_xy.plot([Cx, x0], [Cy, y0], ls="--", color="black", lw=1.2)
ax_xy.scatter(x0, y0, s=90, color="black", zorder=10)

ax_xy.text(Cx+8, Cy+16, r"$(C_x,C_y)$")
ax_xy.text((Cx+x0)/2+5, (Cy+y0)/2, r"$R$")

ax_xy.set_xlabel(r"$x_{\mathrm{TPC}}$ [mm]")
ax_xy.set_ylabel(r"$y_{\mathrm{TPC}}$ [mm]")
ax_xy.set_title(r"Projection: $x_{\mathrm{TPC}}$-$y_{\mathrm{TPC}}$")
ax_xy.set_aspect("equal", adjustable="box")
ax_xy.grid(True, alpha=0.5)

# -----------------------------
# 3D helix
# -----------------------------
ax3d = fig.add_subplot(gs[:, 1], projection="3d")

ax3d.plot(x, y, z, color="black", lw=2.0, alpha=0.8)

for xx, yy, zz, c in zip(xh, yh, zh, colors):
    ax3d.scatter(xx, yy, zz, s=80, color=c, edgecolor="black", linewidth=0.6)

ax3d.scatter(x0, y0, z0, s=120, color="black", zorder=10)

ax3d.set_xlabel(r"$x_{\mathrm{TPC}}$ [mm]", labelpad=8)
ax3d.set_ylabel(r"$y_{\mathrm{TPC}}$ [mm]", labelpad=8)
ax3d.set_zlabel(r"$z_{\mathrm{TPC}}$ [mm]", labelpad=8)
ax3d.set_title(r"Helix track in local coordinates", pad=20)

ax3d.set_xlim(-100, 100)
ax3d.set_ylim(-180, -100)
ax3d.set_zlim(-12, 12)
ax3d.view_init(elev=25, azim=-60)

# -----------------------------
# s-zTPC projection
# -----------------------------
ax_sz = fig.add_subplot(gs[1, 0])

ax_sz.plot(s, z, color="black", lw=1.8)

for ss, zz, c in zip(sh, zh, colors):
    ax_sz.scatter(ss, zz, s=55, color=c, edgecolor="black", linewidth=0.5)

ax_sz.scatter(s0, z0, s=90, color="black", zorder=10)
ax_sz.axhline(Z0, ls="--", color="gray", lw=1)
ax_sz.axvline(0.0, ls="--", color="gray", lw=1)


ax_sz.text(
    55,
    Z0 + 0.4,
    r"$Z_0$",
    fontsize=16
)

ax_sz.text(
    0.07, 0.95,
    r"$z_{\mathrm{TPC}} = Z_0 + D_z s$"
    "\n"
    r"$s = R(\phi-\phi_0)$",
    transform=ax_sz.transAxes,
    ha="left",
    va="top",
    fontsize=16
)

ax_sz.text(
    4,
    z.min()+0.8,
    r"$s=0$",
    fontsize=16
)

ax_sz.set_xlabel(r"$s=R(\phi-\phi_0)$ [mm]")
ax_sz.set_ylabel(r"$z_{\mathrm{TPC}}$ [mm]")
ax_sz.set_title(r"Projection: $s$-$z_{\mathrm{TPC}}$")
ax_sz.grid(True, alpha=0.5)

for ax in [ax_xy, ax_sz]:
    ax.xaxis.set_label_coords(0.92, -0.12)
    ax.yaxis.set_label_coords(-0.12, 0.95)

plt.savefig("helix_coordinate_definition.png", dpi=300, bbox_inches="tight")
#plt.savefig("helix_coordinate_definition.pdf", bbox_inches="tight")
plt.show()

print("s range:", sh.min(), sh.max())
print("zTPC range:", zh.min(), zh.max())
