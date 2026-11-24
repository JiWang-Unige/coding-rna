"""fig_c3rel  (SLOT_WIDE 1664x560)  slide c3-rel

Final-model (T2) gene-chain posteriors on the core-10 panel: (a) pooled reliability (10 equal-width bins),
(b) per-species hit rate of the P >= 0.9 bin against its mean predicted probability.
Source: experiments/C3-CALIB/results/c3_result.json; pooled ECE and per-species ECE checked against report 36 §7.
"""
import json
import os
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

R36 = "docs/review_20261016/36_C3_final_calibration.md"
r = json.load(open(os.path.join(style.REPO, "experiments/C3-CALIB/results/c3_result.json")))
P = r["pooled"]
style.has(R36, style.find(R36, "| **合并** | 53.04"), f"{P['ECE']:.3f}", f"{P['AUROC']:.3f}")
NAME = {"dryas_octopetala": "Dryas", "eschscholzia_californica": "Eschscholzia", "akebia_trifoliata": "Akebia", "hakea_prostrata": "Hakea",
        "zizania_latifolia": "Zizania", "themeda_triandra": "Themeda", "juncus_effusus": "Juncus", "rhynchospora_breviuscula": "Rhynchospora",
        "acorus_gramineus": "Acorus", "amborella_trichopoda": "Amborella"}
for k, v in r["per_species"].items():
    style.has(R36, style.find(R36, f"| {NAME[k]} | "), f"{v['ECE']:.3f}")

W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
a = style.axes_px(fig, 100, 60, 560, 400)
mp = [b["mean_pred"] for b in P["bins"]]
hr = [b["hit_rate"] for b in P["bins"]]
a.plot([0, 1], [0, 1], color=style.GRAY, ls="--", lw=1.6)
a.plot(mp, hr, color=style.BLUE, marker="o", lw=3)
a.set_xlim(0, 1); a.set_ylim(0, 1)
a.set_xlabel("mean predicted P(gene chain correct)")
a.set_ylabel("observed exact-chain hit rate")
a.set_title(f"a  Pooled, 244,268 chains · ECE {P['ECE']:.3f}")
a.text(0.04, 0.97, f"P ≥ 0.9 bin:\npredicted {P['hiconf']['mean_P']:.2f}\nobserved {P['hiconf']['hit_rate']:.2f}",
       ha="left", va="top", fontsize=17, color=style.INK)

b = style.axes_px(fig, 820, 60, 820, 370)
sp = sorted(r["per_species"].items(), key=lambda kv: -kv[1]["hiconf"]["hit_rate"])
xs = range(len(sp))
hits = [v["hiconf"]["hit_rate"] for _, v in sp]
cols = [style.ORANGE if v["hiconf"]["mean_P"] - v["hiconf"]["hit_rate"] > 0.15 else style.BLUE for _, v in sp]
b.bar(xs, hits, color=cols, width=0.66)
for x, (_, v) in zip(xs, sp):
    b.plot([x - 0.38, x + 0.38], [v["hiconf"]["mean_P"]] * 2, color=style.INK, lw=2)
    b.text(x, v["hiconf"]["hit_rate"] - 0.06, f"{v['hiconf']['hit_rate']:.2f}", ha="center", fontsize=15, color=style.BG, fontweight="bold")
b.set_xticks(list(xs), [NAME[k] for k, _ in sp], rotation=35, ha="right", fontsize=16)
b.set_ylim(0.4, 1.0)
b.set_ylabel("hit rate in P ≥ 0.9 bin")
b.set_title("b  P ≥ 0.9 bin by species · bar observed, line predicted")
style.save(fig, "fig_c3rel")
print("ok")
