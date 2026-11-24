"""fig_p4dyn  (SLOT_WIDE 1664x560)  slide c2-dyn

P4 training dynamics: exact-chain F1 of k50U, k50F and T10 (all labels) at the four stops, per seed,
on the 8 new species' validation segments. Source: experiments/C3-CALIB/results/p4_result.json;
each U − F is checked against the report 37 §C.1 table.
"""
import json
import os
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

R37 = "docs/review_20261016/37_P4_dynamics_and_long_gene_failures.md"
res = json.load(open(os.path.join(style.REPO, "experiments/C3-CALIB/results/p4_result.json")))
STEPS = ["768", "1536", "2304", "3072"]
for sd, seedno in (("s1", "1"), ("s2", "2")):
    for st in STEPS:
        r = res["seeds"][sd]["per_step"][st]
        ln = style.find(R37, f"| {seedno} | {st} |")
        style.has(R37, ln, f"{r['diff']:+.2f}".replace("-", "−"), f"{r['U']['F1']:.2f}", f"{r['F']['F1']:.2f}")

W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
for i, (sd, title) in enumerate((("s1", "Seed 1"), ("s2", "Seed 2"))):
    ax = style.axes_px(fig, 110 + i * 800, 90, 680, 380)
    ps = res["seeds"][sd]["per_step"]
    x = list(range(4))
    T = [ps[s]["T10"]["F1"] for s in STEPS]
    U = [ps[s]["U"]["F1"] for s in STEPS]
    F = [ps[s]["F"]["F1"] for s in STEPS]
    ax.plot(x, T, color=style.GRAY, marker="D", ls="--", lw=2)
    ax.plot(x, U, color=style.BLUE, marker="o", lw=3)
    ax.plot(x, F, color=style.ORANGE, marker="o", lw=3, markerfacecolor=style.BG, markeredgewidth=2.6)
    for xi, s in zip(x, STEPS):
        d = ps[s]["diff"]
        ax.text(xi, 82, f"{d:+.1f}", ha="center", fontsize=17, fontweight="bold",
                color=style.BLUE if d > 0 else style.ORANGE)
    ax.set_xticks(x, ["step 768", "1,536", "2,304", "3,072"])
    ax.set_ylim(0, 88)
    ax.set_xlim(-0.4, 3.4)
    ax.set_xlabel("")
    if i == 0:
        ax.set_ylabel("exact-chain F1 (%)")
    ax.set_title(f"{title}   (U − F above each stop)", pad=30)
style.legend_row(fig, 110, 505, [("o", style.BLUE, style.BLUE, "k50U  marginalise"),
                                 ("o", style.BG, style.ORANGE, "k50F  fill background"),
                                 ("dash", style.GRAY, style.GRAY, "T10  all genes labelled (one run)")])
style.save(fig, "fig_p4dyn")
print("ok")
