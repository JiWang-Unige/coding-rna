"""fig_setb  (SLOT_WIDE 1664x560)  slide c1-setb

Claim C1 (zero-shot): trained only on Arabidopsis + rice, PlantCAD2 + LoRA (T2) is 5.4 / 4.9
points below filtered Carbon-A on two time-held-out plants; NT-v2 + LoRA was 24.1 / 29.2 below.

Set B = Dryas octopetala (4 chromosomes, 9,508 reference genes) and Eschscholzia californica
(3 scaffolds, 12,367 genes); neither species is in PlantCAD2 pretraining or in our training
data. Zero-shot, no posterior filter for our models; exact-chain F1 (%); one seed.
Sources (docs/review_20261016/):
  20_M34_stage1_results.md lines 145-146 (T10, T2 with 1-Mb block bootstrap CI), 141 (set
                                         definition); DEV values line 78 (T2) and 121 (T10)
  18_M33_ntv2_followup.md lines 77-78, 81, 88-89, 92 (Carbon-A unfiltered / filtered and
                                         NT-v2 + LoRA, gene-unit bootstrap CI), line 100 (paired gap)
  14_M30_results.md line 79 (Carbon-A unfiltered - filtered, paired)
  15_M31_C1_ntv2_pilot.md lines 118, 124-125 (DEV values of NT-v2 + LoRA and Carbon-A)
"""
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

D = "docs/review_20261016/"
R20, R18, R15 = D + "20_M34_stage1_results.md", D + "18_M33_ntv2_followup.md", D + "15_M31_C1_ntv2_pilot.md"
style.has(R20, 141, "Dryas octopetala 4 条、Eschscholzia californica 3 条", "都不在 PlantCAD2 预训练与 10 个训练物种内")
style.has(R20, 145, "T10（十物种，3072 步）")
style.has(R20, 146, "T2（拟南芥 + 水稻，2304 步）")
t10 = style.cite(R20, 145, "68.11", "66.23", "69.72", "68.50", "67.20", "69.66")
t2 = style.cite(R20, 146, "68.57", "66.55", "70.20", "68.39", "66.98", "69.65")
nt_d = style.cite(R18, 81, "**49.92**", "48.95", "50.95")
nt_e = style.cite(R18, 92, "**44.09**", "43.17", "45.02")
cu_d = style.cite(R18, 77, "72.32", "71.42", "73.19")
cf_d = style.cite(R18, 78, "73.98", "73.08", "74.83")
cu_e = style.cite(R18, 88, "71.75", "70.94", "72.53")
cf_e = style.cite(R18, 89, "73.32", "72.52", "74.10")
style.has(R18, 78, "Carbon-A 过滤")
style.has(R18, 89, "Carbon-A 过滤")
dev = {"T2": style.cite(R20, 78, "**73.88**")[0], "T10": style.cite(R20, 121, "**72.91**")[0],
       "NT": style.cite(R15, 118, "**63.61**")[0], "CU": style.cite(R15, 124, "75.241")[0],
       "CF": style.cite(R15, 125, "77.406")[0]}

rows = [  # label, (dryas, lo, hi), (esch, lo, hi), dev, colour
    ("Carbon-A, filtered", cf_d, cf_e, dev["CF"], style.ORANGE),
    ("Carbon-A, unfiltered", cu_d, cu_e, dev["CU"], style.ORANGE_LIGHT),
    ("T2: PlantCAD2 + LoRA, Arabidopsis + rice", t2[:3], t2[3:], dev["T2"], style.BLUE),
    ("T10: PlantCAD2 + LoRA, ten species", t10[:3], t10[3:], dev["T10"], style.BLUE_LIGHT),
    ("NT-v2 + LoRA, Arabidopsis + rice", nt_d, nt_e, dev["NT"], style.GRAY),
]
# Gaps to filtered Carbon-A. Use the paired values the reports give where they exist (computing
# from rounded F1 values is off by 0.01 in two places); T2 and T10 have point differences only.
R14 = D + "14_M30_results.md"
gaps = [(0.0, 0.0),
        tuple(style.cite(R14, 79, "−1.66", "−1.58")),                 # Carbon-A unfiltered − filtered
        (round(t2[0] - cf_d[0], 2), round(t2[3] - cf_e[0], 2)),        # T2 (point differences)
        (-5.87, -4.82),                                                # T10: report 20 line 157, checked below
        tuple(style.cite(R18, 100, "−24.05", "−29.23"))]               # NT-v2 + LoRA, paired
style.has(R14, 79, "Carbon-A 未过滤 − Carbon-A 过滤")
style.has(R18, 100, "B 未过滤 − Carbon-A 过滤")
style.has(R20, 157, "Z 到 Carbon-A 过滤的差距只有 5.87 / 4.82 点")
assert gaps[2] == (-5.41, -4.93)
style.close(t10[0] - cf_d[0], -5.87); style.close(t10[3] - cf_e[0], -4.82)

# ------------------------------------------------------------------ figure
W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
ax = style.axes_px(fig, 520, 66, 880, 414)
XMIN, XMAX = 40, 80
for i, (label, d, e, dv, col) in enumerate(rows):
    y = 4 - i
    ax.text(XMIN - 0.8, y, label, ha="right", va="center", fontsize=17,
            fontweight="bold" if label.startswith("T2") else "normal", color=style.INK)
    ax.plot([XMIN, XMAX], [y, y], color=style.GRID, lw=1.0, zorder=0)
    for (x, lo, hi), dy, mk in ((d, 0.17, "o"), (e, -0.17, "s")):
        ax.plot([lo, hi], [y + dy, y + dy], color=col, lw=2.2, zorder=2)
        ax.scatter([x], [y + dy], marker=mk, s=105, color=col, edgecolor=col, zorder=3)
    ax.scatter([dv], [y], marker="D", s=80, facecolor=style.BG, edgecolor=style.MUTED,
               linewidth=2.0, zorder=3)
    ax.plot([min(d[0], e[0]), dv], [y, y], color=style.GRAY_LIGHT, lw=1.6, ls=":", zorder=1)
    gd, ge = gaps[i]
    txt = "reference" if i == 0 else f"{gd:+.2f} / {ge:+.2f}".replace("-", "−")
    ax.text(XMAX + 1.0, y, txt, ha="left", va="center", fontsize=17,
            fontweight="bold" if label.startswith("T2") else "normal",
            color=style.BLUE if label.startswith("T2") else style.BODY)
ax.text(XMAX + 1.0, 4.62, "vs filtered Carbon-A\nDryas / Eschscholzia", ha="left", va="center",
        fontsize=15, color=style.MUTED, linespacing=1.1)
ax.set_xlim(XMIN, XMAX)
ax.set_ylim(-0.6, 4.95)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xticks(range(40, 81, 5))
ax.set_xlabel("Exact-chain F1 (%), zero-shot on set B, one seed; lines = 95% CI", labelpad=8)
style.legend_row(fig, 520, 4, [
    ("o", style.MUTED, style.MUTED, "Dryas octopetala"),
    ("s", style.MUTED, style.MUTED, "Eschscholzia californica"),
    ("D", style.BG, style.MUTED, "same model on DEV (upper bound)"),
], gap=34, fontsize=16)
print(style.save(fig, "fig_setb"))
