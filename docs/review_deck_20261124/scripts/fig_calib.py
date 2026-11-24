"""fig_calib  (SLOT_WIDE 1664x560)  slide c3-calib

Claim C3: the partition-function posteriors of NT-v2 + LoRA (M31 B) are calibrated at tau = 1
without rescaling (tau chosen on internal validation), and filtering by posterior gains
+1.7 on DEV and +2.3 to +2.5 on the two set B species.

Raw local results (M33, copied from ont): experiments/M33-NTV2-FOLLOWUP/results/
  selection.json  taus[*].ECE (internal validation, tau grid), selection.tau = 1.0
  dev_report.json calibration["1"].bins / ECE / AUROC (DEV), bootstrap.B_unfiltered.diff_Bfiltered_minus_this
  results_<species>.json posterior_on_B.bins / ECE / AUROC, bootstrap_diffs["NTv2B_post - NTv2B"]
Checked against docs/review_20261016/18_M33_ntv2_followup.md lines 12, 41, 99, 111-113.
Chain-level posterior vs exact-chain hit; 10 equal-width bins; one seed.
"""
import json
import os
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

RES = os.path.join(style.REPO, "experiments/M33-NTV2-FOLLOWUP/results")
R18 = "docs/review_20261016/18_M33_ntv2_followup.md"
sel = json.load(open(os.path.join(RES, "selection.json")))
dev = json.load(open(os.path.join(RES, "dev_report.json")))
dry = json.load(open(os.path.join(RES, "results_dryas_octopetala.json")))
esc = json.load(open(os.path.join(RES, "results_eschscholzia_californica.json")))

assert sel["selection"]["tau"] == 1.0
taus = sorted(sel["taus"], key=float)
ece_tau = [sel["taus"][t]["ECE"] for t in taus]
style.has(R18, 41, "**1.0** | **0.031** | **0.820**", "**−1.400**")
cal = {"DEV": dev["calibration"]["1"], "Dryas": dry["posterior_on_B"], "Eschscholzia": esc["posterior_on_B"]}
style.has(R18, 111, "0.817", "0.031")
style.has(R18, 112, "0.785", "0.043")
style.has(R18, 113, "0.093")
style.close(cal["DEV"]["ECE"], 0.031, 0.0006)
style.close(cal["Dryas"]["ECE"], 0.043, 0.0006)
style.close(cal["Eschscholzia"]["ECE"], 0.093, 0.0006)
g_dev = dev["bootstrap"]["B_unfiltered"]["diff_Bfiltered_minus_this"]
gains = [("DEV", g_dev["point"], *g_dev["CI95"])]
for name, r in (("Dryas", dry), ("Eschscholzia", esc)):
    d = r["bootstrap_diffs"]["NTv2B_post - NTv2B"]
    gains.append((name, d["point"], d["lo"], d["hi"]))
style.has(R18, 12, "+1.68", "+1.31 至 +2.09")
style.has(R18, 99, "+2.47 [+2.12, +2.80]", "+2.31 [+1.95, +2.70]")
for (n, p, lo, hi), ref in zip(gains, [(1.68, 1.31, 2.09), (2.47, 2.12, 2.80), (2.31, 1.95, 2.70)]):
    style.close(p, ref[0]); style.close(lo, ref[1]); style.close(hi, ref[2])

COL = {"DEV": style.BLUE, "Dryas": style.ORANGE, "Eschscholzia": style.GRAY}
MK = {"DEV": "o", "Dryas": "s", "Eschscholzia": "D"}

# ------------------------------------------------------------------ figure
W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()

# (a) reliability
axa = style.axes_px(fig, 92, 62, 410, 400)
axa.plot([0, 1], [0, 1], color=style.MUTED, lw=1.4, ls="--", zorder=1)
for name, c in cal.items():
    xs = [b["mean_pred"] for b in c["bins"]]
    ys = [b["hit_rate"] for b in c["bins"]]
    ns = [b["n"] for b in c["bins"]]
    axa.plot(xs, ys, color=COL[name], lw=2.2, zorder=2)
    axa.scatter(xs, ys, s=[18 + 140 * n / max(ns) for n in ns], marker=MK[name], color=COL[name],
                zorder=3)
axa.set_xlim(0, 1); axa.set_ylim(0, 1)
axa.set_xticks([0, 0.5, 1]); axa.set_yticks([0, 0.5, 1])
axa.set_xlabel("Gene posterior (bin mean)", labelpad=6)
axa.set_ylabel("Fraction exact", labelpad=6)
axa.set_title("Reliability at τ = 1", pad=10)
for i, name in enumerate(cal):
    c = cal[name]
    axa.text(0.97, 0.30 - 0.1 * i, f"{name}: ECE {c['ECE']:.3f}", ha="right", va="center",
             fontsize=15.5, color=COL[name], transform=axa.transAxes)

# (b) tau selection on internal validation
axb = style.axes_px(fig, 640, 62, 380, 400)
tv = [float(t) for t in taus]
axb.plot(tv, ece_tau, color=style.BLUE, lw=2.4, marker="o", markersize=8)
axb.scatter([1.0], [sel["taus"]["1"]["ECE"]], s=220, facecolor="none", edgecolor=style.INK,
            linewidth=2, zorder=4)
axb.annotate("chosen: τ = 1\nECE 0.031", xy=(1.0, sel["taus"]["1"]["ECE"]), xytext=(1.12, 0.36),
             fontsize=16, color=style.INK, ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=style.MUTED, lw=1.3, shrinkB=9))
axb.set_xlim(0, 2.15); axb.set_ylim(0, 0.75)
axb.set_xticks([0.2, 0.5, 1, 1.5, 2])
axb.set_xticklabels(["0.2", "0.5", "1", "1.5", "2"])
axb.set_yticks([0, 0.25, 0.5, 0.75])
axb.set_xlabel("Temperature τ", labelpad=6)
axb.set_ylabel("ECE, internal validation", labelpad=6)
axb.set_title("τ chosen before DEV", pad=10)

# (c) filter gains
axc = style.axes_px(fig, 1290, 62, 340, 400)
for i, (name, p, lo, hi) in enumerate(gains):
    y = 2 - i
    axc.plot([lo, hi], [y, y], color=COL[name], lw=3)
    axc.scatter([p], [y], marker=MK[name], s=120, color=COL[name], zorder=3)
    axc.text(hi + 0.08, y, f"+{p:.2f}", ha="left", va="center", fontsize=17, color=style.INK)
    axc.text(-0.12, y, name, ha="right", va="center", fontsize=17, color=style.INK)
axc.axvline(0, color=style.INK, lw=1.2)
axc.set_xlim(0, 3.4); axc.set_ylim(-0.6, 2.6)
axc.set_yticks([])
axc.spines["left"].set_visible(False)
axc.set_xticks([0, 1, 2, 3])
axc.set_xlabel("Filtered − unfiltered F1, 95% CI", labelpad=6)
axc.set_title("Gain from filtering", pad=10)
print(style.save(fig, "fig_calib"))
