"""fig_panel  (SLOT_WIDE 1664x560)  slide c1-panel

C1 verdict on the independent panel: per-species exact-chain F1 (core 10 plants) for T2 and the baselines,
species sorted by T2 − Carbon-A filtered. Source: experiments/PANEL-P6/results/p6_panel_table.tsv;
medians checked against report 31 §6 (line with '| **中位数** |').
"""
import csv
import os
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

R31 = "docs/review_20261016/31_P6_panel_baselines.md"
rows = list(csv.DictReader(open(os.path.join(style.REPO, "experiments/PANEL-P6/results/p6_panel_table.tsv")), delimiter="\t"))
med = next(r for r in rows if r["species"] == "panel_median")
ln = style.find(R31, "| **中位数** |")
style.has(R31, ln, med["T2"], med["Helixer"], med["Tiberius"], med["ANNEVO"], med["CarbonA_filt"], med["Vipsania_ft"])
sp = [r for r in rows if not r["species"].startswith("panel_")]
assert len(sp) == 10
sp.sort(key=lambda r: float(r["T2"]) - float(r["CarbonA_filt"]))
REFSEQ = {"dryas_octopetala", "eschscholzia_californica", "akebia_trifoliata", "hakea_prostrata", "amborella_trichopoda"}
TOOLS = [("Helixer", style.GRAY_LIGHT, "s"), ("Tiberius", style.GRAY, "s"), ("Vipsania_ft", style.SLATE, "s"),
         ("ANNEVO", style.ORANGE_LIGHT, "s"), ("CarbonA_filt", style.ORANGE, "D"), ("T2", style.BLUE, "o")]

W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
ax = style.axes_px(fig, 330, 70, 1000, 410)
data = sp + [med]
ys = list(range(len(sp) + 1, 1, -1)) + [0]
for y, r in zip(ys, data):
    vals = [float(r[t]) for t, _, _ in TOOLS]
    ax.plot([min(vals), max(vals)], [y, y], color=style.GRID, lw=4, zorder=1)
    for (t, c, m), v in zip(TOOLS, vals):
        ax.scatter([v], [y], marker=m, s=150 if t == "T2" else 80, color=c, zorder=3 if t == "T2" else 2,
                   edgecolor=style.BG if t != "T2" else c, linewidth=0.8)
    name = "Median (10)" if r["species"] == "panel_median" else r["species"].replace("_", " ").capitalize()
    tag = "" if r["species"] == "panel_median" else ("  RefSeq" if r["species"] in REFSEQ else "  GenBank")
    ax.text(24.5, y, name, ha="right", va="center", fontsize=18,
            fontstyle="normal" if tag == "" else "italic", fontweight="bold" if tag == "" else "normal")
    ax.text(81, y, tag.strip(), ha="left", va="center", fontsize=16, color=style.MUTED)
    d = float(r["T2"]) - float(r["CarbonA_filt"])
    ax.text(91, y, f"{d:+.1f}", ha="left", va="center", fontsize=18, color=style.BLUE,
            fontweight="bold" if tag == "" else "normal")
ax.text(81, len(sp) + 1.9, "annotation     T2 − CA filt", ha="left", fontsize=16, color=style.MUTED)
ax.axhline(1.0, color=style.GRID, lw=1.2)
ax.set_xlim(25, 80)
ax.set_ylim(-0.7, len(sp) + 2.2)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xticks([30, 40, 50, 60, 70, 80])
ax.set_xlabel("exact-chain F1 (%), whole panel sequences, one run per tool")
style.legend_row(fig, 330, 0, [(m, c, c, {"CarbonA_filt": "Carbon-A filtered", "Vipsania_ft": "Vipsania ft", "T2": "T2 (ours)"}.get(t, t))
                                 for t, c, m in TOOLS])
style.save(fig, "fig_panel")
print("ok")
