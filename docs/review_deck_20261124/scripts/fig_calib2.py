"""fig_calib2  (SLOT_WIDE 1664x560)  slide c3-pc2

Claim C3 (status): calibration has been measured on four model families and three yardsticks;
the first PlantCAD2 numbers (T10 on set B) are close to NT-v2 on the same species, and
self-training makes the posteriors over-confident. The final model T2 has not been measured.

Sources (docs/review_20261016/):
  09_phase0_phase1_diagnostics.md line 5 (TRAIN = the six training chromosomes), line 32 (our GENERanno posterior AUROC 0.805, E/M pool),
    line 34 (Carbon-A confidence AUROC 0.833), line 69 (ECE 0.156 ours, 0.067 Carbon-A)
  12_M29_A2_temperature.md line 43 (gate model rescaled to tau = 0.03: AUROC 0.834, ECE 0.086)
  18_M33_ntv2_followup.md lines 111-113 (NT-v2 + LoRA at tau = 1: DEV, Dryas, Eschscholzia)
  27_M36_results.md line 53 (T10 on set B; after self-training S1)
Chain-level posterior (or confidence) against exact-chain hit, 10 equal-width bins, one seed.
"""
import sys

from matplotlib import transforms

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

D = "docs/review_20261016/"
gen_auc, = style.cite(D + "09_phase0_phase1_diagnostics.md", 32, "0.805")
ca_auc, = style.cite(D + "09_phase0_phase1_diagnostics.md", 34, "0.833")
gen_ece, ca_ece = style.cite(D + "09_phase0_phase1_diagnostics.md", 69, "ECE 0.156", "ECE 0.067")
gate_auc, gate_ece = style.cite(D + "12_M29_A2_temperature.md", 43, "**0.834**", "**0.086**")
style.has(D + "12_M29_A2_temperature.md", 43, "| **0.03** |")
nt = {}
for ln, key in ((111, "DEV"), (112, "Dryas"), (113, "Esch")):
    c = style.cells(D + "18_M33_ntv2_followup.md", ln)
    nt[key] = (style._num(c[3]), style._num(c[1]))      # ECE, AUROC
assert nt["DEV"] == (0.031, 0.817) and nt["Dryas"] == (0.043, 0.785) and nt["Esch"] == (0.093, 0.760)
style.has(D + "27_M36_results.md", 53, "ECE 0.059 / 0.089、AUROC 0.773 / 0.728",
          "S1 的 ECE 升到 0.187 / 0.208")

rows = [  # label, region, ECE, AUROC, colour, bold
    ("GENERanno posterior, τ = 1 (09)", "TRAIN pool", gen_ece, gen_auc, style.SLATE, False),
    ("GENERanno gate, rescaled τ = 0.03 (M29-A2)", "internal validation", gate_ece, gate_auc, style.SLATE, False),
    ("Carbon-A confidence network (09)", "TRAIN pool", ca_ece, ca_auc, style.ORANGE, False),
    ("NT-v2 + LoRA, τ = 1 (M33)", "DEV", *nt["DEV"], style.BLUE_LIGHT, False),
    ("NT-v2 + LoRA, τ = 1 (M33)", "set B Dryas", *nt["Dryas"], style.BLUE_LIGHT, False),
    ("NT-v2 + LoRA, τ = 1 (M33)", "set B Eschscholzia", *nt["Esch"], style.BLUE_LIGHT, False),
    ("PlantCAD2 T10, τ = 1 (M36 Z)", "set B Dryas", 0.059, 0.773, style.BLUE, True),
    ("PlantCAD2 T10, τ = 1 (M36 Z)", "set B Eschscholzia", 0.089, 0.728, style.BLUE, True),
    ("T10 after self-training (M36 S1)", "set B Dryas", 0.187, None, style.GRAY, False),
    ("T10 after self-training (M36 S1)", "set B Eschscholzia", 0.208, None, style.GRAY, False),
]

# ------------------------------------------------------------------ figure
W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
TOP, HT = 58, 430
axe = style.axes_px(fig, 760, TOP, 470, HT)
axa = style.axes_px(fig, 1300, TOP, 330, HT)
blend = transforms.blended_transform_factory(axe.transAxes, axe.transData)
ys = []
y = len(rows) + 1.0
for i in range(len(rows)):
    if i in (3, 6, 8):
        y -= 0.5
    ys.append(y)
    y -= 1
for y, (label, region, ece, auc, col, bold) in zip(ys, rows):
    axe.text(-(760 - 548) / 470, y, label, ha="right", va="center", fontsize=16,
             fontweight="bold" if bold else "normal", color=style.INK, transform=blend)
    axe.text(-(760 - 562) / 470, y, region, ha="left", va="center", fontsize=15.5,
             color=style.MUTED, transform=blend)
    axe.plot([0, ece], [y, y], color=style.GRID, lw=1.2, zorder=0)
    axe.scatter([ece], [y], s=110, color=col, zorder=3)
    axe.text(ece + 0.006, y, f"{ece:.3f}", ha="left", va="center", fontsize=15,
             color=style.BODY)
    if auc is not None:
        axa.plot([0.5, auc], [y, y], color=style.GRID, lw=1.2, zorder=0)
        axa.scatter([auc], [y], s=110, color=col, zorder=3)
        axa.text(auc + 0.008, y, f"{auc:.3f}", ha="left", va="center", fontsize=15, color=style.BODY)
for ax, xlim, ticks, lab, title in (
        (axe, (0, 0.30), [0, 0.1, 0.2, 0.3], "ECE (lower = better calibrated)", "Calibration"),
        (axa, (0.5, 0.92), [0.5, 0.6, 0.7, 0.8, 0.9], "AUROC (ranking)", "Ranking")):
    ax.set_xlim(*xlim)
    ax.set_ylim(ys[-1] - 0.7, ys[0] + 0.7)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks(ticks)
    ax.set_xlabel(lab, labelpad=6)
    ax.set_title(title, pad=8)
print(style.save(fig, "fig_calib2"))
