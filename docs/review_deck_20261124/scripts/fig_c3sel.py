"""fig_c3sel  (SLOT_WIDE 1664x560)  slide c3-sel

Selective annotation with T2 posteriors on the core-10 panel: (a) pooled precision / recall / F1 when keeping the
top x% of chains by posterior, with Carbon-A's own confidence at 70% and 50% coverage; (b) per-species F1 change
from the frozen threshold log P >= -1.400 (chosen on NT-v2 in M33, not re-tuned for T2).
Source: experiments/C3-CALIB/results/c3_result.json; checked against report 36 §7.
"""
import json
import os
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

R36 = "docs/review_20261016/36_C3_final_calibration.md"
r = json.load(open(os.path.join(style.REPO, "experiments/C3-CALIB/results/c3_result.json")))
P = r["pooled"]
cov = P["coverage"]
style.has(R36, style.find(R36, "| 精度 | 53.4"), *[f"{cov[k]['precision']:.1f}" for k in ("100", "90", "80", "70", "50")])
style.has(R36, style.find(R36, "| **合并** | 53.04"), f"+{P['frozen_thr']['dF1']:.2f}")
NAME = {"dryas_octopetala": "Dryas", "eschscholzia_californica": "Eschscholzia", "akebia_trifoliata": "Akebia", "hakea_prostrata": "Hakea",
        "zizania_latifolia": "Zizania", "themeda_triandra": "Themeda", "juncus_effusus": "Juncus", "rhynchospora_breviuscula": "Rhynchospora",
        "acorus_gramineus": "Acorus", "amborella_trichopoda": "Amborella"}

W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
a = style.axes_px(fig, 100, 60, 640, 400)
ks = ["100", "90", "80", "70", "50"]
x = [int(k) for k in ks]
for key, col, lab in (("precision", style.BLUE, "precision"), ("recall", style.BLUE_LIGHT, "recall"), ("F1", style.INK, "F1")):
    a.plot(x, [cov[k][key] for k in ks], color=col, marker="o", lw=3, label=lab)
ca = P["carbonA"]
a.scatter([70, 50], [ca["precision_at_70"], ca["precision_at_50"]], marker="D", s=90, color=style.ORANGE, zorder=4, label="Carbon-A precision (own conf.)")
a.set_xlim(102, 46)
a.set_ylim(30, 82)
a.set_xticks(x, [f"{k}%" for k in ks])
a.set_xlabel("chains kept, ranked by posterior")
a.set_ylabel("%")
a.set_title("a  Risk–coverage, 10 species pooled")
a.legend(loc="lower left", fontsize=15, ncol=2)

b = style.axes_px(fig, 860, 60, 780, 370)
sp = sorted(r["per_species"].items(), key=lambda kv: -kv[1]["frozen_thr"]["dF1"])
xs = range(len(sp))
d = [v["frozen_thr"]["dF1"] for _, v in sp]
assert min(d) > 0
b.bar(xs, d, color=style.BLUE, width=0.66)
for xi, (_, v) in zip(xs, sp):
    b.text(xi, v["frozen_thr"]["dF1"] + 0.15, f"+{v['frozen_thr']['dF1']:.1f}", ha="center", fontsize=15, color=style.INK)
b.set_xticks(list(xs), [NAME[k] for k, _ in sp], rotation=35, ha="right", fontsize=16)
b.set_ylim(0, 9)
b.set_ylabel("F1 change (points)")
b.set_title(f"b  Frozen threshold: +{P['frozen_thr']['dF1']:.2f} pooled, 10 of 10 species up")
style.save(fig, "fig_c3sel")
print("ok")
