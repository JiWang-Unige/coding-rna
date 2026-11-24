"""fig_dev  (SLOT_WIDE 1664x560)  slide c1-dev

Claim C1 (backbone-agnostic): the same head, grammar and sampling lift every backbone; LoRA
adds ~12 points to both self-supervised backbones; the final PlantCAD2 recipe (T2) reaches 73.88.

DEV = one reused chromosome each of Arabidopsis and rice (in PlantCAD2 pretraining and in the
Carbon-A training data: an upper bound, not a test). Exact-chain F1 (%), pooled over the two
species; 95% CI = 1-Mb block bootstrap where reported; one seed per model.
Sources (all docs/review_20261016/):
  15_M31_C1_ntv2_pilot.md lines 118-125  (NT-v2 frozen / + LoRA, Helixer, Tiberius, ANNEVO,
                                         GENERanno fresh-joint, Carbon-A unfiltered / filtered)
  17_M32_plantcad2_pilot.md line 86     (PlantCAD2 frozen)
  21_M35_plantcad2_lora_result.md line 9 (PlantCAD2 + LoRA, M35)
  20_M34_stage1_results.md line 78      (T2, final recipe)
  18_M33_ntv2_followup.md line 54       (NT-v2 + LoRA + posterior filter)
"""
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

D = "docs/review_20261016/"
ntv2_lora, = style.cite(D + "15_M31_C1_ntv2_pilot.md", 118, "**63.61**")
ntv2_lora_ci = style.cite(D + "15_M31_C1_ntv2_pilot.md", 118, "61.17", "65.62")
ntv2_frozen, = style.cite(D + "15_M31_C1_ntv2_pilot.md", 119, "50.92")
helixer, = style.cite(D + "15_M31_C1_ntv2_pilot.md", 120, "62.851")
tiberius, = style.cite(D + "15_M31_C1_ntv2_pilot.md", 121, "71.805")
annevo, = style.cite(D + "15_M31_C1_ntv2_pilot.md", 122, "75.703")
gener, gener_lo, gener_hi = style.cite(D + "15_M31_C1_ntv2_pilot.md", 123, "70.107", "67.74", "72.11")
ca_unf, = style.cite(D + "15_M31_C1_ntv2_pilot.md", 124, "75.241")
ca_filt, = style.cite(D + "15_M31_C1_ntv2_pilot.md", 125, "77.406")
pc2_frozen, pc2_lo, pc2_hi = style.cite(D + "17_M32_plantcad2_pilot.md", 86, "**59.875**", "56.48", "62.45")
m35, m35_lo, m35_hi = style.cite(D + "21_M35_plantcad2_lora_result.md", 9, "**72.219**", "70.05", "73.97")
t2, t2_lo, t2_hi = style.cite(D + "20_M34_stage1_results.md", 78, "**73.88**", "71.80", "75.63")
filt, filt_lo, filt_hi = style.cite(D + "18_M33_ntv2_followup.md", 54, "**65.29**", "62.95", "67.17")
style.has(D + "15_M31_C1_ntv2_pilot.md", 120, "Helixer 0.3.7")
style.has(D + "15_M31_C1_ntv2_pilot.md", 121, "Tiberius 2.0.8")
style.has(D + "15_M31_C1_ntv2_pilot.md", 122, "ANNEVO 2.3.3")
style.close(ntv2_lora - ntv2_frozen, 12.69, 0.006)    # "LoRA +12.7"
style.close(m35 - pc2_frozen, 12.344, 0.006)          # "+12.34" (report 21 line 11)

# ------------------------------------------------------------------ figure
W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
ax = style.axes_px(fig, 520, 62, 1000, 418)
XMIN, XMAX = 46, 80

B, BL, SL, GR, OR = style.BLUE, style.BLUE_LIGHT, style.SLATE, style.GRAY, style.ORANGE
rows = [  # label, x, lo, hi, colour, filled, marker, bold
    ("PlantCAD2 frozen (M32)", pc2_frozen, pc2_lo, pc2_hi, B, False, "o", False),
    ("PlantCAD2 + LoRA (M35)", m35, m35_lo, m35_hi, B, True, "o", False),
    ("PlantCAD2 + LoRA, final recipe T2", t2, t2_lo, t2_hi, B, True, "*", True),
    ("NT-v2 frozen (M31 A)", ntv2_frozen, None, None, BL, False, "o", False),
    ("NT-v2 + LoRA (M31 B)", ntv2_lora, *ntv2_lora_ci, BL, True, "o", False),
    ("NT-v2 + LoRA + posterior filter (M33)", filt, filt_lo, filt_hi, BL, True, "s", False),
    ("GENERanno preview, frozen (M28, 16 Oct)", gener, gener_lo, gener_hi, SL, False, "o", False),
    ("Helixer 0.3.7", helixer, None, None, GR, True, "D", False),
    ("Tiberius 2.0.8", tiberius, None, None, GR, True, "D", False),
    ("ANNEVO 2.3.3", annevo, None, None, GR, True, "D", False),
    ("Carbon-A, unfiltered", ca_unf, None, None, OR, True, "D", False),
    ("Carbon-A, filtered", ca_filt, None, None, OR, True, "D", False),
]
ys, y = [], 12.6
for i in range(len(rows)):
    if i in (3, 6, 7):
        y -= 0.45
    ys.append(y)
    y -= 1
for y, (label, x, lo, hi, col, filled, mk, bold) in zip(ys, rows):
    ax.text(XMIN - 0.6, y, label, ha="right", va="center", fontsize=16.5,
            fontweight="bold" if bold else "normal",
            color=style.INK if not label.startswith(("Helixer", "Tiberius", "ANNEVO", "Carbon"))
            else style.BODY)
    ax.plot([XMIN, x], [y, y], color=style.GRID, lw=1.0, zorder=0)
    if lo is not None:
        ax.plot([lo, hi], [y, y], color=col, lw=2.4, zorder=2)
    ax.scatter([x], [y], marker=mk, s=330 if mk == "*" else (110 if mk != "D" else 70),
               facecolor=col if filled else style.BG, edgecolor=col, linewidth=2.4, zorder=3)
    xt = (hi if hi is not None else x) + 0.6
    ax.text(xt, y, f"{x:.2f}", ha="left", va="center", fontsize=16.5,
            fontweight="bold" if bold else "normal", color=B if bold else style.BODY)
for ysep in (ys[2] - 0.72, ys[5] - 0.72, ys[6] - 0.72):
    ax.axhline(ysep, color=style.GRID, lw=1.0, xmin=0, xmax=1)
ax.text(XMAX, ys[7] + 0.62, "other tools, own pipelines", ha="right", va="center",
        fontsize=15, color=style.MUTED)
ax.text(XMAX, ys[0] + 0.62, "same head, grammar and sampling", ha="right", va="center",
        fontsize=15, color=style.MUTED)
ax.set_xlim(XMIN, XMAX)
ax.set_ylim(ys[-1] - 0.6, ys[0] + 0.9)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xticks(range(50, 81, 5))
ax.set_xlabel("DEV exact-chain F1 (%), Arabidopsis + rice pooled, one seed; lines = 95% CI where reported",
              labelpad=8)
print(style.save(fig, "fig_dev"))
