"""fig_long  (SLOT_WIDE 1664x560)  slide c4-long

Panel core-10 (247,922 references): (a) exact-chain recall by CDS span, T2 vs Carbon-A unfiltered;
(b) error type of 'gap genes' (T2 misses, Carbon-A hits) for the three largest-gap species and two small-gap controls.
Source: experiments/C3-CALIB/results/long_gene_result.json, checked against report 37 §C.2.
"""
import json
import os
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

R37 = "docs/review_20261016/37_P4_dynamics_and_long_gene_failures.md"
r = json.load(open(os.path.join(style.REPO, "experiments/C3-CALIB/results/long_gene_result.json")))
P = r["pooled"]["by_span"]
LAB = {"<2kb": "< 2 kb", "2-5kb": "2–5 kb", "5-10kb": "5–10 kb", "10-20kb": "10–20 kb", "20-50kb": "20–50 kb", ">=50kb": "≥ 50 kb"}
for k, v in P.items():
    style.has(R37, style.find(R37, f"| {LAB[k]} | "), f"{v['n']:,}", f"{v['recall_CA']:.1f}")

W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
a = style.axes_px(fig, 100, 60, 680, 380)
ks = list(LAB)
xs = range(len(ks))
a.bar([x - 0.2 for x in xs], [P[k]["recall_T2"] for k in ks], width=0.4, color=style.BLUE, label="T2")
a.bar([x + 0.2 for x in xs], [P[k]["recall_CA"] for k in ks], width=0.4, color=style.ORANGE, label="Carbon-A unfiltered")
for x, k in zip(xs, ks):
    a.text(x - 0.2, P[k]["recall_T2"] + 1, f"{P[k]['recall_T2']:.0f}", ha="center", fontsize=15, color=style.BLUE)
a.set_xticks(list(xs), [f"{LAB[k]}\nn={P[k]['n']/1000:.0f}k" if P[k]['n'] >= 10000 else f"{LAB[k]}\nn={P[k]['n']:,}" for k in ks], fontsize=15)
a.set_ylim(0, 80)
a.set_ylabel("exact-chain recall (%)")
a.set_title("a  Recall by CDS span (10 species pooled)")
a.legend(loc="upper right", fontsize=15)

b = style.axes_px(fig, 960, 60, 600, 380)
SPEC = [("hakea_prostrata", "Hakea"), ("amborella_trichopoda", "Amborella"), ("akebia_trifoliata", "Akebia"),
        ("eschscholzia_californica", "Eschscholzia"), ("dryas_octopetala", "Dryas")]
CATS = [("split", "split", style.BLUE), ("single_junction_shift", "one junction shifted", style.BLUE_LIGHT),
        ("complex_chain", "complex", style.SLATE), ("exon_missed", "exon missed", style.GRAY),
        ("absent_no_call", "no call", style.ORANGE), ("other", "other", style.GRAY_LIGHT)]
for y, (k, name) in enumerate(SPEC[::-1]):
    g = r["species"][k]["T2_miss_CA_hit"]
    n = g["n"]
    left = 0
    for c, _, col in CATS:
        v = (sum(val for kk, val in g["categories"].items() if kk not in [x[0] for x in CATS[:-1]]) if c == "other" else g["categories"].get(c, 0)) / n * 100
        b.barh(y, v, left=left, color=col, height=0.62)
        if v >= 9:
            b.text(left + v / 2, y, f"{v:.0f}", ha="center", va="center", fontsize=15, color=style.BG if col not in (style.BLUE_LIGHT, style.GRAY_LIGHT) else style.INK)
        left += v
    b.text(102, y, f"n={n:,}", va="center", fontsize=15, color=style.MUTED)
b.set_yticks(range(len(SPEC)), [n for _, n in SPEC[::-1]], fontsize=17)
b.set_xlim(0, 100)
b.set_xlabel("% of gap genes (T2 misses, Carbon-A hits)")
b.set_title("b  Error type of gap genes")
style.legend_row(fig, 100, 520, [("s", col, col, lab) for _, lab, col in CATS], gap=26, fontsize=15)
style.save(fig, "fig_long")
print("ok")
