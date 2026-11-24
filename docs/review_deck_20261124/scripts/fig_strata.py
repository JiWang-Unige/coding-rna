"""fig_strata  (SLOT_WIDE 1664x560)  slide c4-strata

Claim C4 (descriptive): the remaining gap is in genes with many CDS segments; the backbone
moves the >=6-segment F1 most, the number of training species does not.

DEV exact-chain F1 (%) by number of CDS segments of the reference gene (1 / 2-5 / >=6), one seed.
Sources (docs/review_20261016/):
  15_M31_C1_ntv2_pilot.md lines 118, 119, 123, 124, 125  (NT-v2 B, NT-v2 A, GENERanno, Carbon-A)
  17_M32_plantcad2_pilot.md line 86   (PlantCAD2 frozen)
  21_M35_plantcad2_lora_result.md line 14  (M35)
  20_M34_stage1_results.md line 78 (T2), line 121 (T10)
"""
import re
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

D = "docs/review_20261016/"
TRIPLE = r"(\d+\.\d+) / (\d+\.\d+) / (\d+\.\d+)"


def strata(path, ln, anchor=None, which=0):
    """The strata triple on path:ln: the first one after `anchor`, or the `which`-th one."""
    text = style.line(D + path, ln)
    if anchor:
        text = text[text.index(anchor):]
    found = re.findall(TRIPLE, text)
    assert found, f"{path}:{ln} has no strata triple"
    return [float(x) for x in found[which]]


S = {
    "NT-v2 frozen": strata("15_M31_C1_ntv2_pilot.md", 119),
    "NT-v2 + LoRA": strata("15_M31_C1_ntv2_pilot.md", 118),
    "PlantCAD2 frozen": strata("17_M32_plantcad2_pilot.md", 86, which=-1),   # last triple; the first is CDS-base P/R/F1
    "GENERanno (16 Oct)": strata("15_M31_C1_ntv2_pilot.md", 123),
    "PlantCAD2 + LoRA (M35)": strata("21_M35_plantcad2_lora_result.md", 14),
    "T10 (ten species)": strata("20_M34_stage1_results.md", 121, "1 / 2–5 / ≥6 段"),
    "T2 (final)": strata("20_M34_stage1_results.md", 78, "1 / 2–5 / ≥6 段"),
    "Carbon-A unfiltered": strata("15_M31_C1_ntv2_pilot.md", 124),
    "Carbon-A filtered": strata("15_M31_C1_ntv2_pilot.md", 125),
}
assert S["PlantCAD2 frozen"] == [75.86, 57.79, 50.87]
assert S["T2 (final)"] == [78.56, 72.92, 71.69]
assert S["T10 (ten species)"] == [78.90, 71.87, 69.80]
assert S["NT-v2 + LoRA"] == [75.71, 63.71, 53.98]
assert S["Carbon-A filtered"] == [80.49, 76.23, 76.49]
style.close(S["PlantCAD2 + LoRA (M35)"][2] - S["NT-v2 + LoRA"][2], 14.76)
style.close(S["T2 (final)"][2] - S["NT-v2 + LoRA"][2], 17.71)

# ------------------------------------------------------------------ figure
W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()

# (a) slope across strata
axa = style.axes_px(fig, 90, 64, 520, 410)
LINES = [("Carbon-A filtered", style.ORANGE, "-", 2.6), ("T2 (final)", style.BLUE, "-", 3.2),
         ("T10 (ten species)", style.BLUE_LIGHT, "--", 2.4), ("GENERanno (16 Oct)", style.SLATE, "-", 2.2),
         ("NT-v2 + LoRA", style.GRAY, "-", 2.4)]
xs = [0, 1, 2]
for name, col, ls, lw in LINES:
    ys = S[name]
    axa.plot(xs, ys, color=col, ls=ls, lw=lw, marker="o", markersize=8, zorder=3)
    axa.text(2.08, ys[2], f"{name}  {ys[2]:.1f}", ha="left", va="center", fontsize=15.5,
             color=col if col not in (style.GRAY, style.BLUE_LIGHT) else style.BODY,
             fontweight="bold" if name.startswith("T2") else "normal")
axa.set_xticks(xs)
axa.set_xticklabels(["1", "2–5", "≥6"])
axa.set_xlim(-0.15, 2.0)
axa.set_ylim(50, 82)
axa.set_xlabel("CDS segments in the reference gene", labelpad=6)
axa.set_ylabel("DEV exact-chain F1 (%)", labelpad=6)
axa.set_title("The gap opens with exon count", pad=10)

# (b) >=6-segment ladder
axb = style.axes_px(fig, 1240, 64, 390, 410)
order = ["NT-v2 frozen", "PlantCAD2 frozen", "NT-v2 + LoRA", "GENERanno (16 Oct)",
         "PlantCAD2 + LoRA (M35)", "T10 (ten species)", "T2 (final)", "Carbon-A unfiltered",
         "Carbon-A filtered"]
COL = {"NT-v2 frozen": style.GRAY_LIGHT, "NT-v2 + LoRA": style.GRAY, "PlantCAD2 frozen": style.BLUE_LIGHT,
       "GENERanno (16 Oct)": style.SLATE, "PlantCAD2 + LoRA (M35)": style.BLUE_LIGHT,
       "T10 (ten species)": style.BLUE_LIGHT, "T2 (final)": style.BLUE,
       "Carbon-A unfiltered": style.ORANGE_LIGHT, "Carbon-A filtered": style.ORANGE}
for i, name in enumerate(order):
    y = len(order) - 1 - i
    v = S[name][2]
    axb.barh(y, v, height=0.62, color=COL[name], zorder=2)
    axb.text(v + 0.8, y, f"{v:.1f}", ha="left", va="center", fontsize=15.5, color=style.INK,
             fontweight="bold" if name.startswith("T2") else "normal")
    axb.text(-1.5, y, name, ha="right", va="center", fontsize=15.5, color=style.INK,
             fontweight="bold" if name.startswith("T2") else "normal")
axb.set_xlim(0, 85)
axb.set_ylim(-0.6, len(order) - 0.4)
axb.set_yticks([])
axb.spines["left"].set_visible(False)
axb.set_xticks([0, 20, 40, 60, 80])
axb.set_xlabel("DEV F1, genes with ≥6 segments", labelpad=6)
axb.set_title("≥6 segments, by backbone", pad=10)
print(style.save(fig, "fig_strata"))
