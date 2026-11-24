"""fig_fail  (SLOT_WIDE 1664x560)  slide c4-replay

Claim C4 (mechanism, descriptive): errors are wrong joins more than missed genes, and most
missed genes are representable and scored positive but lose to a competing chain.

(a) Set B failure modes of every reference gene (exact chain / overlapping but not exact /
    no overlap), NT-v2 + LoRA (M31 B, M33 evaluation) vs filtered Carbon-A, Dryas and
    Eschscholzia pooled (gene-weighted), by number of CDS segments.
    Raw local file: experiments/M33-NTV2-FOLLOWUP/results/diag_failure_modes.json
    Report: docs/review_20261016/18_M33_ntv2_followup.md line 14 (33.6% / 35.0% overall).
(b) Matched final-update score replay of an M28 GENERanno model with joint encoding (joint24,
    DEV 66.623, report 04 line 191; not the final fresh-joint model, 70.107) on TRAIN:
    docs/review_20261016/04_M28_current_state.md line 131 (7,762 misses: 2,486 negative energy,
    5,067 positive but out-competed, 209 not supported by the grammar); primary report
    experiments/M28/full_grammar/reports/full_grammar_matched_energy_diagnostic.md.
"""
import json
import os
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

FM = json.load(open(os.path.join(style.REPO, "experiments/M33-NTV2-FOLLOWUP/results/diag_failure_modes.json")))
R04 = "docs/review_20261016/04_M28_current_state.md"
R18 = "docs/review_20261016/18_M33_ntv2_followup.md"
style.has(R18, 14, "33.6% / 35.0%")
for sp, ref in (("dryas_octopetala", 0.336), ("eschscholzia_californica", 0.350)):
    style.close(FM[sp]["NTv2B"]["all"]["overlap_not_exact"], ref, 0.0006)
neg, pos, unsup = style.cite(R04, 131, "7,762", "2,486", "5,067", "209")[1:]
style.has(R04, 131, "7,762 条漏检中 2,486 条能量为负、5,067 条为正但被竞争链压过、209 条语法不支持")
assert neg + pos + unsup == 7762


def pooled(model, stratum):
    out, ntot = {"exact": 0.0, "overlap_not_exact": 0.0, "missed": 0.0}, 0
    for sp in ("dryas_octopetala", "eschscholzia_californica"):
        r = FM[sp][model][stratum]
        ntot += r["n"]
        for k in out:
            out[k] += r[k] * r["n"]
    return {k: v / ntot for k, v in out.items()}, ntot


rows = []
for stratum, lab in (("1", "1 segment"), ("2-5", "2–5 segments"), (">=6", "≥6 segments")):
    for model, mlab in (("CarbonA_filt", "Carbon-A filtered"), ("NTv2B", "NT-v2 + LoRA")):
        fr, n = pooled(model, stratum)
        rows.append((lab, mlab, fr, n))
ov6 = {r[1]: r[2]["overlap_not_exact"] for r in rows if r[0] == "≥6 segments"}
style.close(ov6["NT-v2 + LoRA"], 0.5817, 0.0006)
style.close(ov6["Carbon-A filtered"], 0.1636, 0.0006)

# ------------------------------------------------------------------ figure
W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
CAT = [("exact", "exact chain", style.BLUE), ("overlap_not_exact", "overlaps, not exact", style.ORANGE_LIGHT),
       ("missed", "no overlap (missed)", style.GRAY_LIGHT)]

axa = style.axes_px(fig, 230, 104, 700, 372)
ys = [7.0, 6.2, 4.3, 3.5, 1.6, 0.8]
for y, (lab, mlab, fr, n) in zip(ys, rows):
    left = 0
    for key, _, col in CAT:
        w = fr[key] * 100
        axa.barh(y, w, left=left, height=0.62, color=col, zorder=2)
        if w >= 7:
            axa.text(left + w / 2, y, f"{w:.0f}%", ha="center", va="center", fontsize=15.5,
                     color="#FAF9F5" if key == "exact" else style.INK)
        left += w
    axa.text(-2, y, mlab, ha="right", va="center", fontsize=16, color=style.INK)
for ytop, (lab, n) in zip((7.0, 4.3, 1.6), [(r[0], r[3]) for r in rows[::2]]):
    axa.text(0, ytop + 0.52, f"{lab} · {n:,} reference genes", ha="left", va="center",
             fontsize=15.5, color=style.MUTED)
axa.set_xlim(0, 100)
axa.set_ylim(0.3, 7.75)
axa.set_yticks([])
axa.spines["left"].set_visible(False)
axa.set_xticks([0, 25, 50, 75, 100])
axa.set_xlabel("Reference genes on set B (Dryas + Eschscholzia), %", labelpad=6)


axb = style.axes_px(fig, 1150, 104, 480, 372)
parts = [("positive energy, lost to a competing chain", pos, style.BLUE),
         ("negative energy (scored below background)", neg, style.GRAY),
         ("not allowed by the grammar", unsup, style.GRAY_LIGHT)]
for i, (lab, v, col) in enumerate(parts):
    y = 2 - i
    axb.barh(y, v, height=0.42, color=col, zorder=2)
    axb.text(v + 120, y, f"{int(v):,} ({v / 7762:.0%})", ha="left", va="center", fontsize=16.5,
             color=style.INK)
    axb.text(0, y + 0.47, lab.replace("\n", " "), ha="left", va="center", fontsize=15.5,
             color=style.BODY)
axb.set_xlim(0, 8200)
axb.set_ylim(-0.5, 2.75)
axb.set_yticks([])
axb.spines["left"].set_visible(False)
axb.set_xticks([0, 2000, 4000, 6000])
axb.set_xticklabels(["0", "2,000", "4,000", "6,000"])
axb.set_xlabel("TRAIN reference genes missed (7,762)", labelpad=6)
axb.set_title("Why TRAIN genes are missed (M28)", pad=34, fontsize=17)

axa.set_title("Set B: where each reference gene ends up", pad=34, fontsize=17)
style.legend_row(fig, 230, 2, [("s", c, c, t) for _, t, c in CAT], gap=34, fontsize=16)
print(style.save(fig, "fig_fail"))
