"""Build the 24 November 2026 review deck from the 2 November one.

usage: python build_1124.py <deck project dir>      (e.g. <scratch>/deck/project, a copy of review_deck_20261102/deck_source/project)

- patches stale sentences in the slides that are kept (exact-string replacements; each must match once),
- writes the new and rewritten slides,
- rewrites deck.json (order, sections, title).
Every number below was checked against the repository (verification TSVs in review_deck_20261124/verification/).
Page numbers are placeholders; run the skill's renumber.py afterwards.
"""
import json
import os
import sys

P = sys.argv[1]
S = os.path.join(P, "slides")

BG, PANEL, INK, BODY, MUTED, BLUE, RUST, LINE = "#FAF9F5", "#F2F0EA", "#1B2430", "#3D4757", "#6B7380", "#1F5F8B", "#A14E1C", "#D9D6CC"
SEC = f"background:{BG};color:{INK};font-family:'IBM Plex Sans', Arial, sans-serif;padding:104px 128px 160px;display:flex;flex-direction:column;gap:28px"


def head(kicker, title):
    return (f'<p style="font-size:24px;font-weight:600;letter-spacing:3px;text-transform:uppercase;color:{BLUE}">{kicker}</p>\n'
            f'<h2 style="font-family:\'Source Serif 4\', Georgia, serif;font-size:46px;font-weight:600;line-height:1.18;color:{INK};'
            f'padding-bottom:22px;border-bottom:2px solid {LINE}">{title}</h2>\n')


def foot(text):
    return (f'<div style="position:absolute;left:128px;bottom:64px;width:1664px;display:flex;flex-direction:row;justify-content:space-between">\n'
            f'<p style="font-size:24px;color:{MUTED}">{text}</p>\n<p style="font-size:24px;color:{MUTED}">0</p>\n</div>\n')


def write(sid, body, notes, bg=BG):
    assert len(notes) <= 4000, (sid, len(notes))
    html = f'<section id="{sid}" data-transition="fade" style="{SEC.replace(BG, bg, 1)}">\n{body}<aside>{notes}</aside>\n</section>\n'
    open(os.path.join(S, f"{sid}.html"), "w").write(html)


def figure_slide(sid, kicker, title, blob, alt, footer, notes):
    write(sid, head(kicker, title) + f'<img src="/_blob/{blob}" alt="{alt}" style="width:1664px;height:560px;object-fit:contain">\n' + foot(footer), notes)


def table(header, rows, widths, size=24, num_cols=()):
    th = "".join(f'<th style="padding:8px 14px;width:{w}%;text-align:{"right" if i in num_cols else "left"}">{h}</th>' for i, (h, w) in enumerate(zip(header, widths)))
    out = f'<table style="font-size:{size}px;color:{INK};font-variant-numeric:tabular-nums">\n<tr style="background:{PANEL}">{th}</tr>\n'
    for r in rows:
        out += "<tr>" + "".join(f'<td style="padding:7px 14px;text-align:{"right" if i in num_cols else "left"}">{c}</td>' for i, c in enumerate(r)) + "</tr>\n"
    return out + "</table>\n"


BADGE = {"Done": f"color:{BG};background:{BLUE};padding:2px 16px", "Paused": f"color:{BG};background:{RUST};padding:2px 16px",
         "Needs your decision": f"color:{RUST};background:{BG};border:2px solid {RUST};padding:0 14px",
         "Queued": f"color:{BLUE};background:{BG};border:2px solid {BLUE};padding:0 14px",
         "Not started": f"color:{MUTED};background:{BG};border:2px solid {MUTED};padding:0 14px"}


def cards(items, badge=True, size=24):
    out = '<div style="display:grid;grid-template-columns:1fr 1fr 1fr;gap:20px">\n'
    for it in items:
        out += f'<div style="display:flex;flex-direction:column;gap:8px;background:{BG if badge else PANEL};border:2px solid {LINE};border-radius:8px;padding:20px 24px">\n'
        if badge:
            out += f'<p style="align-self:start;font-size:24px;font-weight:600;{BADGE[it[0]]};border-radius:999px">{it[0]}</p>\n'
            it = it[1:]
        out += f'<h3 style="font-size:28px;font-weight:600;color:{INK}">{it[0]}</h3>\n'
        if len(it) == 3:
            out += f'<p style="font-size:24px;line-height:1.35;color:{MUTED}">{it[1]}</p>\n'
        out += f'<p style="font-size:{size}px;line-height:1.35;color:{BODY}">{it[-1]}</p>\n</div>\n'
    return out + "</div>\n"


def note_box(text, color=BLUE):
    return f'<p style="font-size:26px;line-height:1.4;color:{BODY};border-left:6px solid {color};padding:6px 0 6px 24px">{text}</p>\n'


def patch(sid, pairs, new_id=None):
    p = os.path.join(S, f"{sid}.html")
    h = open(p).read()
    for a, b in pairs:
        assert h.count(a) == 1, (sid, a[:80], h.count(a))
        h = h.replace(a, b)
    if new_id:
        h = h.replace(f'id="{sid}"', f'id="{new_id}"', 1)
        os.remove(p)
        p = os.path.join(S, f"{new_id}.html")
    open(p, "w").write(h)


# ===================================================================================================== cover
write("cover", f'''<div style="position:absolute;left:0;top:0;width:24px;height:1080px;background:{BLUE}"></div>
<p style="font-size:24px;font-weight:600;letter-spacing:3px;text-transform:uppercase;color:{BLUE}">Systematic review · 24 November 2026</p>
<h1 style="font-family:'Source Serif 4', Georgia, serif;font-size:76px;font-weight:600;line-height:1.12;color:{INK};max-width:1500px">A gene grammar for DNA language models: what held, what failed</h1>
<p style="font-family:'Source Serif 4', Georgia, serif;font-size:36px;font-style:italic;line-height:1.3;color:{BODY};padding-bottom:36px;border-bottom:2px solid {LINE}">What we did since June, the problems we hit, the results that did not work, and the claims that survive</p>
<div style="display:flex;flex-direction:column;gap:8px">
<p style="font-size:32px;font-weight:600;color:{INK}">Ji Wang</p>
<p style="font-size:26px;color:{MUTED}">University of Geneva · coding-rna project</p>
</div>
<p style="position:absolute;right:128px;bottom:64px;width:1100px;text-align:right;font-size:24px;color:{MUTED}">Draft v2 · repository main, 24 November · updates the 2 November review</p>
''', """Purpose. This deck is for you, not for an audience: it lays out the whole project so that we can decide together what story the paper tells. It updates the review of 2 November (28 slides). Since then the decisive tests came in: the ten-plant panel (C1), a second seed for P4 (C2), an outside review of the claims, a test of C2 on a real incomplete annotation (M42), and the final analyses of the posteriors (C3) and of long genes (C4).

How to read it. The blue kicker says which part a slide belongs to. The title states the conclusion with its key number. The footer gives the source (report and section), the evaluation region and the number of runs. Every number was checked against the repository on 24 November by four read-only verification passes (verification/ folder); figures are redrawn from result files with asserts on the report lines.

Order. Slides 2–6: the problem, how we got here, the problems we hit, what changed, and the claims map. Then the evidence for each claim (C1–C4), the list of what was tested and dropped, the integrity checks, and finally status, pending branches, limits and the questions I want to settle with you.

Yardsticks. F1 is exact-chain F1: a predicted gene counts only if every CDS segment matches the reference. DEV = one reused development chromosome each of Arabidopsis and rice (an upper bound: both are in PlantCAD2 pretraining and Carbon-A training). Set B = held-out chromosomes of Dryas octopetala and Eschscholzia californica, scored once per model. Panel = ten independent plant genomes, scored once per tool on 17 November.""")

# ===================================================================================================== journey
J = [
    ("Jun–Sep", "M1–M27: per-base classifiers, fixed callers, structure heads, three backbones", "Most failed cross-species bars; the structure head was a no-go"),
    ("29 Sep–14 Oct", "M28: full CDS grammar CRF, exact whole-contig decoding (GENERanno)", "DEV 70.11, below Tiberius 71.81 and ANNEVO 75.70"),
    ("16–17 Oct", "M29: residual head on frozen Carbon-A outputs", "Below Carbon-A; accuracy race dropped"),
    ("17–22 Oct", "M30: 25% labels on Carbon-A features · M31–M35: NT-v2, PlantCAD2, LoRA", "U − F +16.5 / +14.2; PlantCAD2 + LoRA 72.22 chosen"),
    ("22 Oct–2 Nov", "M34: 2 vs 10 species, half-labelled arms · M36: self-training on the target genome", "More species: no gain · M36 NO-GO · P4 holds +3.24"),
    ("3–17 Nov", "Independent ten-plant panel, six baseline tools", "C1 fails: median 55.05 vs 67.65 for filtered Carbon-A"),
    ("17–23 Nov", "M41 intron-length states · P4 second seed · outside review", "M41 no effect · P4 holds (+9.10) · two critiques of C2"),
    ("23–24 Nov", "M42 oracle-free C2 · C3 calibration · P4 dynamics · long genes", "M42 stops · posteriors rank but over-confident · long genes split"),
]
write("journey", head("How we got here", "Five months, four pivots: from trying to beat Carbon-A to a grammar layer whose value is in how it learns, not in its accuracy")
      + table(["When", "What we tried", "Outcome"], [(f"<b>{a}</b>", b, c) for a, b, c in J], [14, 50, 36], size=24)
      + foot("Reports 03–37 · git log · verification/D_history.md · numbers are exact-chain F1 unless stated"),
      """The arc in one paragraph. From June to September (M1–M27) we tried per-base classifiers with fixed callers, structure heads and several backbones; most failed the cross-species bars we had set. In late September the full CDS grammar (M28) gave a working exact-chain annotator, but on DEV it was below Tiberius and ANNEVO (70.11 against 71.81 and 75.70). On 16–17 October a residual head on frozen Carbon-A outputs (M29) stayed below Carbon-A, and we stopped trying to win on accuracy. The story moved to what the grammar layer can do that the others cannot: train with unknown positions marginalised (M30), and give gene posteriors (M33).

Backbones. NT-v2 + LoRA reached 63.61 on DEV but collapsed zero-shot on set B (49.92 / 44.09). PlantCAD2 + LoRA (M35, 72.22; final recipe T2 73.88) closed most of that: set B 68.57 / 68.39. Adding eight training species (T10) did not help (−0.96 DEV, −0.14 set B). Self-training on the target genome (M36) made every arm worse (−1.5 to −7.7).

The two decisive tests of November. The ten-plant panel (17 November) failed the C1 criterion: T2's median 55.05 against 67.65 for filtered Carbon-A. P4, the test of C2, held in both seeds (+3.24, +14.97). Then an outside review (ChatGPT, 23 November, checked by me) pointed out that the C2 tests placed the unknown mask using the complete annotation. The oracle-free test M42 (24 November) stopped at gate 1: on a real incomplete annotation, partial-label training is 2.0 points below the zero-shot model.

What survived. A grammar layer that trains on any of three backbones; marginalising unknown positions is safer than filling them when the missing positions are known; posteriors that rank gene chains well (AUROC 0.77) though over-confident; and a clear failure mode (genes longer than 10 kb are split). The following slides give each of these with its evidence.""")

# ===================================================================================================== problems
write("problems", head("How we got here · problems", "Most delays were compute: one shared GPU with TE priority, a frozen cluster, a fragile host, and two of our own bugs")
      + cards([
          ("One shared GPU", "ont's RTX 4090 is shared with the TE project, which has priority. Our queue was paused on 27 Oct, 7 Nov and since 23 Nov (M39)."),
          ("Frozen cluster", "Baobab credits fell below threshold in October: no jobs until 1 December. Modal stayed with TE. All baselines ran on ont."),
          ("Fragile host", "ont restarted on 4, 5 and 18 Nov; CPU cores 0 and 8 segfault; the disk filled twice; OOMs in Carbon-A and Vipsania."),
          ("Our bugs", "A gate read an empty TE-claim file as busy: 17 h idle (fixed 21 Nov). Missed waiter notices left the GPU idle 2.7 h."),
          ("New cluster (CAU)", "Prebuilt mamba kernels need GLIBC 2.32, CAU has 2.28: built from source. M42 gate 1 cost 41.5 yuan."),
          ("Design surprises", "Arm L collapsed after 8 steps; P4's margin depends on the checkpoint; DEV was never a test of generalisation."),
      ], badge=False)
      + foot("Reports 00, 20, 29, 31, 33, 35 · commits 7f28143, f185c89, 0a4bb60, 01e56f6 · memory notes · verification/D_history.md"),
      """Compute (the main constraint). One RTX 4090 on ont, shared with the TE project, which has priority by your rule. Pauses: M34 against M36 on 27 October (SIGSTOP/SIGCONT); TE claimed the GPU from 7 November for about 2.5 days (the panel gate was changed to yield); since 23 November TE's timing benchmark has the GPU and M39 (insects) is stopped at step 2,304 of I10_k50U, expected to resume between 29 November and 2 December. Other users' jobs also blocked k50F (5.4 GB free against a 7.5 GB gate) and k25U.

Baobab. Credits were 147,263 on 16 October, under the 200k threshold, so no Baobab jobs until the billing year resets on 1 December (asked first each time). Modal is reserved for TE. All six panel baselines therefore ran on ont (about 85 GPU hours) between 9 and 15 November.

Host problems. ont crashed on 4 November (21:55) and shut down on 5 November (05:17), and rebooted hard on 18 November during k50F_s2 step 4. Cores 0 and 8 of the i9-13900K segfault, so they were removed from our CPU sets (11–12 November). The disk filled on 10 and 13–14 November (other users), pausing panel fine-tuning. Carbon-A crashed and retried; Rhynchospora ran out of memory; Vipsania needed batch size 1.

Our own bugs. The M41 gate tested the TE claim with -e (exists) instead of -s (non-empty), so an empty claim file kept the queue waiting 17 hours (fixed in 7f28143; the panel gate had the same bug, fixed on 14 November). Three waiter notifications were missed on 9–10 November, leaving the GPU idle 1.7, 0.5 and 0.5 hours. On 24 November two CAU jobs wrote the same schedule file at once (fixed).

CAU (24 November). You allocated the CAU HPC account hhbgroup02 to coding-rna. The prebuilt mamba-ssm 2.3.0 and causal-conv1d 1.6.0 wheels need GLIBC 2.32 and CAU has 2.28, so both were compiled from source with the same versions; numba and pandas were matched to ont's environment. M42 gate 1 cost about 41.5 yuan (16.85 GPU hours plus about 1.1 yuan of smoke tests and builds).

Design surprises. M42's arm L (per-position cross-entropy with the CRF decoder) collapsed to zero predicted genes after 8 training steps: a small per-base shift is multiplied by gene length in the CRF score. We found this before any evaluation and dropped the arm (your decision). P4's margin depends heavily on the checkpoint (C2 slides). DEV is in both PlantCAD2's pretraining and Carbon-A's training, so it only ranks models.""", bg=BG)

# ===================================================================================================== changes
write("changes", head("What changed since 2 November", "C1 failed on the panel, C2 survived only where the missing genes are known, and C3 and C4 have their final numbers")
      + table(["", "2 November", "24 November"], [
          ("<b>C1</b> zero-shot", "waiting for the panel", "fails: median 55.05 vs 67.65 (−12.60)"),
          ("<b>C2</b> P4", "+3.24, one seed", "holds in two seeds (+3.24, +14.97); unstable"),
          ("<b>C2</b> in reality", "not tested", "M42 stops: U − zero-shot −2.01"),
          ("<b>C3</b> posteriors", "NT-v2 only", "T2 panel: AUROC 0.77, ECE 0.088, over-confident"),
          ("<b>C4</b> errors", "segment strata on DEV", "genes ≥ 10 kb: recall 0–24% vs 46–59%"),
          ("Seeds", "one run per model", "T2 three runs; P4 two runs per arm"),
          ("Target", "NC if C1 holds", "Genome Biology (your choice, 24 Nov)"),
      ], [20, 30, 50])
      + foot("Reports 20 §16, 31 §6, 35 §8, 36 §7, 37 · presubmission plan 20261123 · dashboard decisions of 24 Nov"),
      """Seven changes, in the order they arrived.

1. C1 (17 November, report 31 §6). The pre-registered panel criterion needed T2's median on the core-10 panel to be at least the medians of Helixer, Tiberius and Vipsania fine-tuned, and within 6.0 points of filtered Carbon-A. Result: T2 55.05 [54.23, 55.78]; Helixer 42.10, Tiberius 59.25, Vipsania ft 57.41, ANNEVO 63.08, Carbon-A unfiltered 65.73 and filtered 67.65 [67.05, 68.25]; T2 − filtered Carbon-A = −12.60 [−13.65, −11.65]. Fails. The Nature Communications route closed; the outline of 17 November moved C2 to the lead.

2. P4 seed 2 (21 November, report 20 §16.2). Registered on 9 November before any seed-2 output: the two-seed mean of k50U − k50F ≥ +2.0 and the per-species mean positive in ≥ 6 of 8. Seed 2: 69.03 − 54.06 = +14.97, 8 of 8. Mean +9.10, 8 of 8. Holds. Most of seed 2's margin is k50F_s2 under-calling (recall 43.7%).

3. Outside review (23 November). ChatGPT read the repository; I verified two critiques: the unknown mask in M30/P4 was placed with the complete annotation (report 19 §3.10 already called it an upper bound), and T2 trained with zero labels of the eight new species scores 68.84 on them, against k50U 67.86 / 69.03.

4. M42 (24 November). Route B, gates 0 and 1, your decision. Gate 0 passed; gate 1 stopped: U − Z = −2.01 [−2.56, −1.50]. You then chose route A plus a merge analysis.

5. C3 (report 36). T2's posteriors on the panel: AUROC 0.768, ECE 0.088; the P ≥ 0.9 bin is right 77% of the time.

6. C4 (report 37). Recall of genes with CDS span 10–20, 20–50 and ≥ 50 kb: T2 24.1, 3.3, 0.0; Carbon-A 59.2, 53.6, 45.6.

7. Seeds and target. T2 now has three runs (DEV 73.88 / 74.09 / 71.65); P4 two per arm. You chose Genome Biology on 24 November and asked for this review before deciding the story.""")

# ===================================================================================================== claims
write("claims", head("Claims map", "Each claim survives in a narrower form; the accuracy claim (C1 zero-shot) and the real-annotation test of C2 failed")
      + table(["Claim", "What survives", "What failed", "Slides"], [
          ("<b>C1</b> One grammar layer, any backbone", "Trains on NT-v2, PlantCAD2 and GENERanno; the backbone sets accuracy (DEV, controlled)", "Zero-shot near Carbon-A: panel −12.60", "C1"),
          ("<b>C2</b> Training on incomplete annotation", "Marginalising beats filling when missing genes are known: M30, P4 in 2 seeds", "Real old → new annotation (M42): no gain over zero-shot", "C2"),
          ("<b>C3</b> Gene-level posteriors", "Rank chains well (AUROC 0.77); frozen threshold helps 10 of 10 species", "Calibration across species: over-confident", "C3"),
          ("<b>C4</b> Where errors are", "Genes ≥ 10 kb are split; errors are wrong joins more than missed genes", "Intron-length states (M41) did not fix it", "C4"),
      ], [24, 38, 28, 10], size=24)
      + foot("Outline 20261117 · presubmission plan 20261123 · reports 20, 31, 35, 36, 37 · status of 24 November"),
      """This slide is the table of contents and the honest scorecard.

C1, as registered on 27 October: "the same grammar layer works on several backbones, and with a plant-pretrained backbone it approaches Carbon-A zero-shot on new plants". The first half holds on DEV under a controlled comparison (same head, grammar, sampling): NT-v2 + LoRA 63.61, PlantCAD2 + LoRA 72.22, GENERanno 70.11; final T2 73.88. The second half failed on the panel (median gap −12.60; also below Tiberius by 4.20 and Vipsania ft by 2.36).

C2: "partial-observation training (marginalise unknown positions) beats filling them with background". Holds in M30 (Carbon-A features, +16.50 / +14.15) and P4 (PlantCAD2, two seeds +3.24 / +14.97, 8 of 8 species each). Both placed the unknown mask using the complete annotation, so they are an upper-bound setting. M42, the realistic version (old RefSeq release as labels, protein evidence as the unknown mask, new release as truth), stopped at gate 1: U is 2.01 points below the zero-shot model and only beats fill-background on missing-locus recall (+2.22). M36 earlier showed the same with pseudo-labels.

C3: "the partition-function model gives calibrated gene posteriors for selective annotation". The ranking part holds on the final model across the panel; the calibration part does not (pooled ECE 0.088 is under 0.10, but five species miss the P ≥ 0.9 bin by more than 0.15, which by the rule written before computing must be stated in the main text).

C4: descriptive. The panel gap to Carbon-A is concentrated in long genes: half of the recall gap comes from the 8.4% of genes with CDS span ≥ 10 kb, mostly through splitting. M41, the one grammar change aimed at this, had no detectable effect.

What would refute each claim now: C1 backbone part, a backbone comparison where the ranking flips with seeds; C2 upper-bound part, M39 insects (paused) failing P4-I; C3 ranking part, another evaluation where the frozen threshold lowers F1; C4, a long-gene analysis on another truth.""")

# ===================================================================================================== C1 kept slides
patch("c1-dev", [
    ("One seed for every model; the +2 seeds of T2 are queued.", "T2 now has three runs: DEV 73.88 / 74.09 / 71.65 (mean 73.21 ± 1.35, report 20 §14); every other model one run."),
    ("(slides 7 and 9 do)", "(set B and the panel slides do)"),
    ("The rule should be applied to every row before the panel is scored.", "The panel scoring (report 31) applied the rule to every tool and also reports the rule-off value."),
])
patch("c1-setb", [
    ("(next slide)", "(ten species did not help either; see Tested and dropped)"),
    ("Two species are not a panel; the claim needs the independent panel (slide 9).", "Two species are not a panel; on the ten-plant panel the gap to filtered Carbon-A widened to 12.60 points (the verdict slide)."),
    ("with “approaches” defined by the panel criterion on slide 24, not by these two species.", "this sentence is withdrawn: the panel criterion failed on 17 November."),
])
patch("c1-panel", [
    ("Ten plants and 4.79 Gb, none in PlantCAD2’s pretraining as a species; half with RefSeq and half with submitter annotation",
     "The panel: ten plants, 4.79 Gb, none in PlantCAD2’s pretraining as a species; half RefSeq, half submitter-annotated"),
    ("· nothing scored yet", "· scored once, 17 Nov"),
    ("Things to keep in mind when the numbers come.", "Things to keep in mind when reading the numbers."),
    ("The rule (slide 24) counts", "The rule (next slide) counts"),
    (">C1 · The test that decides<", ">C1 · The panel<"),
], new_id="c1-panelset")
figure_slide("c1-panel", "C1 · The test that decided",
             "C1 fails: on the ten-plant panel T2’s median is 55.05, 12.60 below filtered Carbon-A and below Tiberius and Vipsania",
             "7d347652b5c3682f85a14d2258105b14",
             "Dot plot of exact-chain F1 for ten panel species and six tools. Medians: Helixer 42.10, T2 55.05, Vipsania ft 57.41, Tiberius 59.25, ANNEVO 63.08, filtered Carbon-A 67.65. T2 minus filtered Carbon-A per species: Amborella −28.0, Hakea −21.7, Akebia −17.1, Zizania −5.8, Dryas −5.2, Eschscholzia −5.0, Juncus −4.3, Acorus −4.2, Rhynchospora −3.5, Themeda −3.1; median −12.6.",
             "Report 31 §6 · core-10 panel, whole sequences · exact-chain F1, CDS ≥ 96 bp · one run per tool, scored once on 17 Nov",
             """The rule (registered 27 October, report 28 §6; tiers in outline §0 and report 31 §2.4). On the core-10 panel, zero-shot, no posterior filter: C1 holds if T2's species median is at least the medians of Helixer a0080, Tiberius and Vipsania fine-tuned (added by your decision vipsania-c1-criterion on 4 November, before any Vipsania number), and within θ = 6.0 points of filtered Carbon-A. Partly holds: gap 6–10, or below one of those tools. Fails: gap over 10.

Result (17 November, scored once). Medians: T2 55.05 [54.23, 55.78]; Helixer 42.10; Tiberius 59.25; Vipsania ft 57.41 (pre-trained 54.69); ANNEVO 63.08; Carbon-A unfiltered 65.73, filtered 67.65 [67.05, 68.25]. T2 − filtered Carbon-A = −12.60 [−13.65, −11.65]; T2 − Tiberius −4.20; T2 − Vipsania ft −2.36. Verdict: fails.

Where the gap is. Three RefSeq species with long, many-exon genes: Amborella −28.0 [−28.96, −27.11], Hakea −21.7, Akebia −17.1. On the five submitter-annotated (GenBank) species the medians are T2 56.61, Tiberius 53.32, Vipsania ft 56.39, filtered Carbon-A 62.45 (gap 5.8); on the five RefSeq species T2 53.49 against 70.69 (gap 17.2). The C4 long-gene slides explain most of the RefSeq gap.

The seed rule (registered 7 November) was triggered: T2_s2 and T2_s3 should be run on the panel. They have not run; they would only report the spread and do not change the verdict.

What it means for the paper. The accuracy claim is gone; zero-shot numbers become an honest benchmark and a limitation. The outline of 17 November moved C2 to the lead and set Genome Biology or NAR Genomics and Bioinformatics as targets.

Caveats. One run per tool. Dryas and Eschscholzia contain the set B chromosomes that T2 had been scored on. Amborella is in Carbon-A's training data (older assembly). GeneCAD (PlantCAD2 + CRF, the closest published method) is missing from the table; whether to run it is on the questions slide.""")

# ===================================================================================================== C2
patch("c2-m30", [("(next slide)", "(the P4 slides)")])
figure_slide("c2-p4", "C2 · Learning from incomplete annotation",
             "P4 holds in both seeds: marginalising beats fill-background in 8 of 8 species, by +3.24 and +14.97 (mean +9.10)",
             "0b2dcc3235fbf9d01cec5102dc122c96",
             "Horizontal bars of k50U minus k50F per species for seed 1 and seed 2, with the two-seed mean. Pooled: +3.24 and +14.97, mean +9.10. Every species positive in both seeds, from Dioscorea +0.20 / +8.20 to Daucus +4.68 / +25.84.",
             "Report 20 §13, §16.2 · 8 new species, validation regions, 5,793 references · each arm at its selected step 3,072 · two runs per arm",
             """The test. Ten training species with PlantCAD2 + LoRA: Arabidopsis and rice fully labelled, and eight new species with half of their genes labelled. k50U treats the other half's positions as unknown and sums over them; k50F calls them background. Scored on the eight new species' internal-validation regions (5,793 reference chains), each arm at the checkpoint its own rule selects on Arabidopsis + rice (step 3,072 in all four arms).

Criteria. Seed 1 (registered before results, plan 19): k50U − k50F ≥ +2.0 pooled, same sign in ≥ 6 of 8. Seed 2 (registered 9 November, before any seed-2 output): the two-seed mean ≥ +2.0 and the per-species mean positive in ≥ 6 of 8.

Results. Seed 1: 67.86 − 64.62 = +3.24, 8 of 8 (2 November). Seed 2: 69.03 − 54.06 = +14.97, 8 of 8 (21 November). Mean +9.10, per-species means +4.20 to +15.26, 8 of 8. Holds.

What drives the margin. In seed 1, recall: k50U finds 4,023 of the 5,793 reference chains, k50F 3,609. In seed 2 almost all of the margin is k50F_s2 under-calling (recall 43.7% against 62.3% in seed 1); k50U_s2 itself is only 1.2 points above k50U. So the average is driven by one checkpoint; the paper should report the two seeds separately.

Against full labels (T10, 69.62, one run): two-seed means k50U −1.18, k50F −10.28. Against no labels at all: T2, which never saw these species, scores 68.84 (three T2 runs: 68.84 / 67.28 / 65.83). This is the second critique on the next slides.

On DEV (fully labelled in both arms) fill-background still costs: k50F − k50U = −1.29 and −0.92, both CIs excluding 0.

Caveat. Both seeds placed the unknown positions using the complete annotation (which genes were unlabelled is known exactly), so this is the upper-bound setting.""")
figure_slide("c2-dyn", "C2 · Learning from incomplete annotation",
             "But training is unstable: U − F is positive at 3 of 4 stops in each seed, negative at step 2,304, and F1 can fall 40 points between stops",
             "d303a32361ccd0e4e97a5f2a158e6506",
             "Two line charts, seed 1 and seed 2, F1 of k50U, k50F and T10 at steps 768, 1536, 2304, 3072. U minus F: seed 1 +39.3, +6.1, −3.9, +3.2; seed 2 +48.3, +6.1, −2.5, +15.0. In seed 2 both arms collapse to about 18–21 at step 2304; T10 rises from 54.0 to 69.6.",
             "Report 37 §C.1 (plan written before computing) · 8 new species, validation regions · CIs from a 1-Mb block bootstrap",
             """What this adds (report 37 part A, descriptive, plan committed before computing). The P4 verdict uses each arm's selected checkpoint. Here are all four stops of both seeds, with the difference split into a "true positives" term (2Δtp over F's pred + ref) and a "number of predictions" term.

Seed 1, U − F [95% CI], tp term / pred term, species positive: 768 +39.33 [+36.88, +41.69], +65.26 / −25.93, 8 of 8; 1,536 +6.12 [+4.33, +7.93], +6.27 / −0.15, 8 of 8; 2,304 −3.94 [−5.02, −2.93], +3.09 / −7.04, 0 of 8; 3,072 +3.24 [+2.29, +4.18], +7.41 / −4.17, 8 of 8.
Seed 2: 768 +48.27 [+45.41, +51.30], 8 of 8; 1,536 +6.11 [+5.15, +7.11], 8 of 8; 2,304 −2.49 [−4.13, −0.99], 1 of 8; 3,072 +14.97 [+13.08, +16.86], 8 of 8.

Reading. (1) Where U wins, it wins on true positives: it finds more correct genes; F's under-calling (negative pred term) actually shrinks the margin. (2) Step 2,304 is negative in both seeds: in seed 1 because U over-calls (7,204 predictions against 5,909 for F); in seed 2 because both arms collapse (F1 18.1 and 20.6). (3) Partial-label training is unstable: within one arm F1 moves by 40–50 points between neighbouring stops (k50F 11–20 at step 768; k50U_s2 18 at 2,304), while full-label T10 rises steadily (53.99, 65.92, 68.92, 69.62). (4) Pairing every stop of U with every stop of F (16 combinations): seed 1 positive in 10, range −13.98 to +48.01; seed 2 positive in 13, range −40.78 to +57.69.

How to write it. "Marginalising let the model find more correct genes at most checkpoints, but partial-label training was unstable" — not a single effect size. The mean +9.10 is a by-product of two checkpoints.

Count check: all 160 region-level counts (5 arms × 4 stops × 8 regions) equal metrics.jsonl.""")
patch("problem", [("Slide 3 shows the full list of changes.", "The changes slide lists what moved since."), ("(slide 13)", "(C2 slides)")])
patch("method", [("(slide 21)", "(integrity · statistics)")])
patch("c2-m36", [("That revision is still waiting for the user’s decision (outline-m36-review, slide 25).", "The outline of 17 November adopted that revision.")])
write("c2-critique", head("C2 · Two critiques (23 November)", "An outside review found two weak points: the unknown mask was placed with the full truth, and partial labels added nothing over zero-shot")
      + '<div style="display:flex;flex-direction:row;gap:28px">\n'
      + f'<div style="flex:1;display:flex;flex-direction:column;gap:14px;background:{PANEL};border-radius:8px;padding:26px 30px">\n<h3 style="font-size:32px;font-weight:600;color:{INK}">1 · The mask borrowed the truth</h3>\n'
      + f'<p style="font-size:26px;line-height:1.4;color:{BODY}">In M30 and P4 the positions of the unlabelled genes came from the complete annotation and were then marked unknown.</p>\n'
      + f'<p style="font-size:26px;line-height:1.4;color:{BODY}">A real genome does not say where its missing genes are. These tests measure an upper bound; report 19 §3.10 already called it that.</p>\n</div>\n'
      + f'<div style="flex:1;display:flex;flex-direction:column;gap:14px;background:{PANEL};border-radius:8px;padding:26px 30px">\n<h3 style="font-size:32px;font-weight:600;color:{INK}">2 · No net gain from partial labels</h3>\n'
      + f'<p style="font-size:26px;line-height:1.4;color:{BODY}">On the eight new species: T2 with <b>no</b> labels of them 68.84; k50U 67.86 / 69.03; full labels (T10) 69.62.</p>\n'
      + f'<p style="font-size:26px;line-height:1.4;color:{BODY}">What C2 shows is that marginalising avoids the loss of filling (k50F 64.62 / 54.06), not that partial labels help.</p>\n</div>\n</div>\n'
      + note_box("<b style=\"color:#1B2430\">Response:</b> M42 tests C2 on a real incomplete annotation (next slides). Caveat on critique 2: T2’s three runs span 65.83–68.84 on these regions, so the comparison depends on which T2 run is used.")
      + foot("Presubmission plan 20261123 §1 · reports 19 §3.10, 20 §13–16 · checked by me against the reports"),
      """Where this came from. On 23 November, at your request, the repository was pushed and ChatGPT ("看论文" project) assessed whether the results make a Nature-family story. I verified every point before writing it down (presubmission plan §1). Two critiques were correct and material.

Critique 1. The U arm in M30 and P4 takes the complete reference annotation, removes a fraction of genes from the labels, and marks exactly those genes' loci as unknown. In a newly sequenced genome nobody knows where the missing genes are: they sit somewhere in the unannotated background. So M30 and P4 test "partial observation when the missing positions are known". Our own report 19 §3.10 had labelled this an upper-bound setting; the paper outline had not.

Critique 2. If partial labels are useful, a model trained with them should beat a model that never saw the species. On the same eight validation regions: T2 (trained only on Arabidopsis + rice, zero labels of the eight species) 68.84; k50U 67.86 (seed 1) and 69.03 (seed 2); T10 (all labels) 69.62; k50F 64.62 and 54.06. So half the labels plus marginalisation is roughly as good as no labels at all, while half the labels plus filling is clearly worse. The claim that survives is defensive: marginalising prevents a loss.

Verification note (24 November): T2's three seeds give 68.84 / 67.28 / 65.83 on these regions (mean 67.32); the two-seed k50U mean is 68.45, 1.1 above that T2 mean. So "no net gain" holds against the best T2 run and becomes a small gain against the average run: the comparison is within seed variation either way.

Response chosen. Route B, gates 0 and 1 (your decision, 23 November): build an oracle-free test where the labels and the unknown mask come only from information available in practice. That is M42 (next three slides). Literature checked at the same time: partial-label CRF training is not new (Tsuboi et al., COLING 2008, the same objective); CRFs for gene prediction exist (Conrad, CRAIG); GeneCAD already combines PlantCAD2 with a CRF. Novelty can only come from what the combination enables on real incomplete annotation.""")
write("c2-m42design", head("C2 · Oracle-free test (M42)", "M42 tests C2 on a real incomplete annotation: the 2016 RefSeq release of jujube as labels, the 2018 release as truth")
      + '<div style="display:flex;flex-direction:row;gap:36px">\n'
      + f'<div style="flex:1;display:flex;flex-direction:column;gap:12px">\n'
      + "".join(f'<p style="font-size:25px;line-height:1.35;color:{BODY}"><b style="color:{INK}">{a}</b> {b}</p>\n' for a, b in [
          ("Data:", "Ziziphus jujuba, one assembly, two releases (AR100 2016 → AR101 2018); 6 chromosomes train (A), 6 evaluate (E)"),
          ("Labels:", "10,472 complete old genes on A as positives"),
          ("Unknown:", "cross-species protein hits ±2 kb (V1), or everything (V3)"),
          ("Arms:", "Z zero-shot T2 · U marginalise · F fill · U3 all unknown; 512 LoRA steps from T2"),
          ("Rules:", "written before results (protocol 35); new release read only at the end")])
      + "</div>\n<div style=\"flex:1\">\n"
      + table(["Gate 0 audit on A", "Missing loci fully unknown", "Unknown share"], [
          ("V1 evidence ±2 kb", "<b>32.8%</b> (247 / 753)", "24.4%"), ("Random, like V1", "11.2%", "21.2%"), ("V2 evidence ±10 kb", "57.1%", "56.5% ✗"),
          ("V3 all unknown", "81.7%", "100%"), ("Oracle (as in P4)", "100%", "1.7%")], [40, 34, 26], size=24)
      + "</div>\n</div>\n"
      + note_box("Gate 0 passed with V1: it covers missing genes 2.9× better than chance, but two thirds of them still sit partly in background.")
      + foot("Protocol 35 §1–§7 · GCF_000826755.1 · audit on training chromosomes A only · 753 missing loci on A"),
      """Design (protocol 35, committed before any data were downloaded). Species chosen from NCBI metadata only: same assembly accession with an older and a newer RefSeq annotation release, chromosome level, not in PlantCAD2's pretraining, enough newly added genes. First choice Ziziphus jujuba (Rosales; AR100 March 2016 → AR101 May 2018; 438 Mb); confirmation species Solanum pennellii; an insect added for gate 2 by your decision.

Split. Twelve chromosomes ranked by length; odd ranks train (A), even ranks evaluate (E).
Positives. Old-release representative transcripts on A that are complete and valid (start/stop, frame, canonical splice sites, no gaps): 10,472 of 11,271 genes.
Evidence. miniprot alignments of the same ten species' proteins used in M36 (not jujube's own): hits with ≥ 50% coverage and ≥ 40% identity.
Masks frozen before the new release was opened: V1 = evidence spans ±2 kb; V2 = ±10 kb; V3 = every non-positive base. Positives ±1 kb are never unknown.

Gate 0 audit (A only, new release read after freezing). 753 genes are new in AR101 (no CDS overlap with AR100); 1,734 revised; 8,996 identical; positives cover 91.7% of identical genes. Missing loci fully inside the unknown region: V1 32.8% (random mask of V1's size 11.2%), V2 57.1% but 56.5% of non-positive bases unknown (rule ≤ 30%: fails), V3 81.7%, oracle 100% with only 1.7% unknown. Rule: V1 passes (≥ 2× random, ≤ 30% unknown) and goes to gate 1; V3 runs as the no-negative-constraint reference.

Gate 1 arms. Z = T2 zero-shot. U = positives + V1 unknown + rest background. F = positives + everything else background. U3 = positives + V3. L (per-position cross-entropy) was dropped after it collapsed in the smoke test, before any scoring (your decision). Continue rule after that revision: U − Z > 0 and U − F > 0 in F1, U's missing-locus recall above Z's, and U's precision not more than 2 points below Z's.

Compute. CAU ChaoBaiHe 4090 (hhbgroup02), about 41.5 yuan in total.""")
figure_slide("c2-m42", "C2 · Oracle-free test (M42)",
             "Gate 1 stops: on real missing genes U is 2.0 points below zero-shot T2 and finds no more of them; it beats fill-background only on missing loci",
             "9926d6180d1a7887e5dc330b77aa53e6",
             "Two bar charts on jujube evaluation chromosomes. F1: Z 59.75, U 57.74, F 58.11, U3 55.28; old annotation used as prediction 74.68. Missing-locus recall (675 loci): Z 16.15, U 15.70, F 13.48, U3 19.70.",
             "Protocol 35 §8.3 · jujube E chromosomes, 10,266 references, 675 missing loci · one run per arm",
             """Result (24 November, scored once on the Mac; the new release never left it). References on E: 10,266 CDS-assessable primaries of AR101, of which 675 are missing from AR100 (274 with protein evidence, 401 without), 7,848 identical and 1,743 revised.

F1 / precision / recall / missing-locus recall (with evidence / without) / recall of old-identical genes / predictions overlapping no gene:
Z 59.75 / 60.10 / 59.41 / 16.15 (22.3 / 12.0) / 72.2 / 407.
U 57.74 / 57.24 / 58.25 / 15.70 (20.1 / 12.7) / 71.0 / 477.
F 58.11 / 60.56 / 55.85 / 13.48 (16.4 / 11.5) / 68.4 / 282.
U3 55.28 / 50.76 / 60.68 / 19.70 (21.9 / 18.2) / 73.5 / 1,501.
Old annotation used as a prediction: 74.68 / 74.46 / 74.90 / 0 / 98.0 / 225.

Paired differences, 1-Mb block bootstrap (154 blocks, 2,000 draws): U − Z F1 −2.01 [−2.56, −1.50], missing-locus recall −0.44 [−1.98, +1.20]. U − F −0.37 [−0.75, +0.01] and +2.22 [+0.91, +3.62]. U − U3 +2.46 [+2.04, +2.95] and −4.00 [−5.75, −2.35].

Verdict: stop. None of the three continue conditions holds (U − Z and U − F both ≤ 0; U's missing-locus recall 15.70 ≤ 16.15; precision 2.86 below Z).

Reading (descriptive). (1) 512 steps of adaptation with the old annotation lowered F1 in all three trained arms (−1.6 to −4.5), as in M36. (2) The only detectable advantage of marginalising over filling is on missing loci (+2.2). (3) Without negatives (U3) the model finds more missing loci (19.7%, even 18.2% without evidence) at the cost of 21% more predictions. (4) Every model is far below simply keeping the old annotation (74.68): in practice a model could only add new loci.

Consequence (protocol 35 §4.4). C2 is written as "partial observation when the missing positions are known" (upper bound); route A; no gate 2.""")
write("c2-merge", head("C2 · Could it help in practice? (post hoc)", "Added on top of the old annotation, model genes in its gaps recover 7–20% of missing loci but lower F1, because only 5–15% of them are right")
      + table(["Old annotation plus", "Added chains", "Right among added", "F1", "Missing loci found"], [
          ("nothing", "—", "—", "<b>74.68</b>", "0"),
          ("Z inside evidence regions", "431", "14.8%", "73.75", "9.5%"),
          ("U inside evidence regions", "465", "12.3%", "73.57", "8.4%"),
          ("F inside evidence regions", "309", "14.9%", "74.01", "6.8%"),
          ("U3 inside evidence regions", "891", "7.6%", "72.21", "10.1%"),
          ("Z anywhere in the gaps", "949", "11.5%", "72.40", "16.2%"),
          ("U anywhere in the gaps", "1,140", "9.3%", "71.73", "15.7%"),
          ("U3 anywhere in the gaps", "2,492", "5.3%", "67.77", "19.7%"),
      ], [34, 16, 18, 14, 18], size=24, num_cols=(1, 2, 3, 4))
      + foot("Protocol 35 §8.4 · rules written before running · jujube E · old-identical genes kept at 98.0% in every row"),
      """Why this analysis (your decision "路线A加合并分析", 24 November). If the model cannot replace an old annotation, can it at least add the genes the old annotation missed, without harming the genes it already has right? The rules were written and committed (e2478e9) before running.

Rules. Keep every old CDS chain on E. Add a model chain if no CDS base overlaps any old CDS on either strand (G0, "anywhere in the gaps"), or additionally only if its CDS span lies entirely inside the E evidence region, i.e. protein hits ±2 kb (G1, "inside evidence regions"; 64.6 Mb, 43% of E). Same ruler and bootstrap as gate 1.

Results. Every merged set keeps 98.0% of the old-identical genes. Missing-locus recall rises to 6.8–19.7%, but the added chains are mostly wrong (5.3–14.9% hit a new-release chain exactly), so F1 drops 0.7 to 6.9 points below the old annotation. The pre-written "practical value" condition (F1 not lower and missing-locus recall higher) is not met.

Paired differences (1-Mb blocks): old + U (G1) − old: F1 −1.11 [−1.29, −0.92], missing loci +8.44 [+5.99, +10.90]. old + U (G1) − old + Z (G1): −0.19 [−0.29, −0.08] and −1.04 [−2.00, −0.15]. old + U (G1) − old + F (G1): −0.44 and +1.63 [+0.78, +2.59]. old + U (G0) − old: −2.94 and +15.70. old + Z (G1) − old: −0.92 and +9.48.

Reading. The evidence gate halves the additions and slightly raises their precision, but what it finds is almost only evidence-supported loci (0.7–2.2% of the evidence-free ones). U is again no better than the zero-shot model; it only beats F on missing loci (+1.63), consistent with gate 1. For the paper: an honest upper limit on recovering missed genes, and a stated limitation that model calls in annotation gaps are of low precision.""")

# ===================================================================================================== C3
patch("c3-calib", [("This is NT-v2, not the final PlantCAD2 model (next slide).", "This is NT-v2, not the final PlantCAD2 model; the next two slides show the final model on the panel.")])
figure_slide("c3-rel", "C3 · Final model on the panel",
             "T2’s posteriors rank chains well but are over-confident everywhere: the P ≥ 0.9 bin is right 77% of the time, and 58% in Themeda",
             "23b5110e4e722d5afb7a452becdbb7ba",
             "Left: pooled reliability curve, observed hit rate below predicted in every bin; P at least 0.9 bin predicted 0.94, observed 0.77; ECE 0.088. Right: observed hit rate of the P at least 0.9 bin per species against predicted about 0.94: Dryas 0.86, Eschscholzia 0.82, Akebia 0.82, Amborella 0.81, Zizania 0.81, Hakea 0.77, Juncus 0.75, Rhynchospora 0.74, Acorus 0.70, Themeda 0.58.",
             "Report 36 §7 (plan written before computing) · core-10 panel, 244,268 chains · τ = 1 · exact-chain hit · one run",
             """What was done (report 36; plan committed 8d2f82d before computing). Every gene chain T2 predicted on the core-10 panel has a posterior P(g) = exp(log_p) at τ = 1 (the temperature M33 chose). A chain is a hit if its exact CDS chain is in the panel reference (same ruler as report 31). Counts were checked first: per-species F1 of T2 and Carbon-A equal the panel table exactly.

Pooled reliability (10 equal-width bins). Observed hit rate is below the mean prediction in every bin, and the gap grows with confidence: 0.45 → 0.42, 0.75 → 0.67, 0.94 → 0.77. Ranking is monotone. ECE 0.088 [0.084, 0.091]; Brier 0.199.

Per species. ECE 0.047 (Dryas) to 0.219 (Themeda) and 0.214 (Acorus). The P ≥ 0.9 bin is predicted at about 0.94 in every species; observed: Dryas 0.86, Eschscholzia 0.82, Akebia 0.82, Amborella 0.81, Zizania 0.81, Hakea 0.77, Juncus 0.75, Rhynchospora 0.74, Acorus 0.70, Themeda 0.58.

The rule written before computing (report 36 §6): if pooled ECE > 0.10, or any species' high-confidence bin is more than 0.15 below its prediction, the main text must say "calibration is not reliable across species". Pooled ECE 0.088 does not trigger it; five species do (exact gaps: Themeda 0.353, Acorus 0.249, Juncus 0.199, Rhynchospora 0.197, Hakea 0.169). So: not reliable across species, systematically over-confident.

History. M33 (NT-v2 + LoRA): DEV ECE 0.031. T10 on set B (M36 §4): 0.059 / 0.089. The final model on new species is in the same range as T10 on set B, worse than NT-v2 on DEV.

Not done. T2 posteriors on DEV (would need GPU). A temperature refit on T2's internal validation (CPU) could reduce ECE but would not fix the species-to-species spread.""")
figure_slide("c3-sel", "C3 · Final model on the panel",
             "Ranking works: a threshold frozen on another backbone raises F1 in 10 of 10 species (+3.08 pooled); AUROC 0.77 against 0.58 for the decode gain",
             "69dbe860b88b3e23b8e5f35d020ed682",
             "Left: risk-coverage, keeping 100, 90, 80, 70, 50 percent of chains by posterior: precision 53.4, 58.8, 64.0, 67.9, 73.1; F1 53.0, 55.3, 56.4, 55.4, 48.3; Carbon-A precision with its own confidence 73.9 at 70 percent and 78.4 at 50 percent. Right: F1 change from the frozen threshold per species, all positive: Hakea +7.7, Amborella +7.5, Akebia +5.0, Acorus +1.9, Zizania +1.7, Dryas +1.7, Rhynchospora +1.5, Eschscholzia +1.3, Juncus +1.3, Themeda +1.0.",
             "Report 36 §7 · threshold log P ≥ −1.400 from M33 (NT-v2), not re-tuned for T2 · Carbon-A ranked by its own confidence",
             """Ranking. AUROC of log P against exact-chain hit: 0.768 pooled (species 0.682–0.874), against 0.581 for the decoder's gain score; difference +0.187 [+0.182, +0.192] (1-Mb block bootstrap, 4,779 blocks). AP 0.736 against 0.572. Carbon-A's dedicated confidence network ranks its own predictions with AUROC 0.751 (ECE 0.079); the two AUROCs are on different prediction sets, so not strictly paired.

Selective annotation. Keeping chains by posterior: precision 53.4 → 58.8 → 64.0 → 67.9 → 73.1 at 100 / 90 / 80 / 70 / 50% kept; F1 peaks at 80% (56.4 against 53.0). At 70% coverage T2's precision is 6.0 points below Carbon-A's at the same coverage [−6.3, −5.7]; this mostly reflects T2's weaker base predictions (unfiltered F1 7.4 points lower), not the posterior.

Frozen threshold. log P ≥ −1.400 was chosen by M33 on NT-v2's internal validation and never re-tuned. On the panel it keeps 85% of T2's chains and raises F1 in 10 of 10 species: Hakea +7.68, Amborella +7.52, Akebia +5.01, then +0.97 to +1.90; pooled 53.04 → 56.12 (+3.08). Carbon-A's own conf ≥ 0.1 filter gains +1.61 (60.41 → 62.02). Filtered T2 (56.12) remains below filtered Carbon-A (62.02).

What C3 can say. The gene-level posterior is a useful ranking score, far better than the decoder's gain and comparable to Carbon-A's separately trained confidence network; a threshold set on another backbone transferred to ten independent species. The probabilities themselves are over-confident and should not be read as error rates (previous slide).""")

# ===================================================================================================== C4
patch("c4-replay", [("(slide 16)", "(previous slide)")])
figure_slide("c4-long", "C4 · Where the errors are",
             "Half of T2’s recall gap to Carbon-A is genes longer than 10 kb: recall falls to 24%, 3% and 0%, and the main error is splitting",
             "06735b9be30453b994200fe0f179ed9e",
             "Left: exact-chain recall by CDS span, T2 versus unfiltered Carbon-A: under 2 kb 49 vs 51, 2–5 kb 65 vs 68, 5–10 kb 53 vs 64, 10–20 kb 24 vs 59, 20–50 kb 3 vs 54, at least 50 kb 0 vs 46. Right: error type of genes T2 misses and Carbon-A finds: split 44 percent in Hakea, 39 in Amborella, 31 in Akebia, 8–9 in Dryas and Eschscholzia, where one shifted junction dominates.",
             "Report 37 §C.2 (plan written before computing) · core-10 panel, 247,922 references · M29 R4 error rules · one run per tool",
             """What was done (report 37 part B). For every panel reference (247,922, same ruler as report 31), whether T2 and unfiltered Carbon-A reproduce its exact chain, stratified by CDS span and by exon count; each miss is put into one structural error class with the M29 R4 rules (no call, wrong strand only, merge, split, endpoint only, exon missed, exon extra, one junction shifted, complex).

Recall. Pooled T2 52.7, Carbon-A 59.1 (gap 6.5). By CDS span, T2 / Carbon-A: < 2 kb 49.1 / 50.9 (106,293 genes); 2–5 kb 65.4 / 67.9 (86,000); 5–10 kb 53.0 / 64.0 (34,891); 10–20 kb 24.1 / 59.2 (13,831); 20–50 kb 3.3 / 53.6 (5,800); ≥ 50 kb 0.0 / 45.6 (1,107). Genes ≥ 10 kb are 8.4% of the references but contribute 3.3 of the 6.5-point recall gap. By exon count the gap grows gradually (1: 50.2 / 52.7; 6–9: 56.9 / 65.8; 20+: 41.1 / 59.7); by span it falls off a cliff near 10 kb.

Gap genes (T2 misses, Carbon-A hits; 25,900 pooled; the reverse 9,794). Pooled: one junction shifted 23%, split 23%, single-exon miss 12%, no call 11%, complex 11%. Hakea (4,493): split 44%; Amborella (4,184): 39%; Akebia (3,850): 31%. Controls Dryas (1,650) and Eschscholzia (2,105): one junction shifted 35% / 32%, split only 8–9%. Monocots Zizania and Themeda: single-exon misses (29%) and no calls.

Reading. In the three species with the largest gap, long genes are cut into pieces: a long-range joining problem, consistent with the older C4 finding (errors are wrong joins). The grammar and decoder have no length cap (exact skyline decoding over whole sequences). A hypothesis, not tested: training windows have a 24,576-bp core, so genes longer than that were never seen complete in training. Short genes (2–5 kb) differ by only 2.5 points, mostly single shifted splice sites.""")
write("c4-m41", head("C4 · A fix that did not work", "Separate long- and short-intron states (M41) did not close the multi-exon gap: −0.74 on genes with ≥ 6 segments, and 12% over-calling on set B")
      + table(["", "T2", "M41", "M41 − T2 [95% CI]"], [
          ("DEV, all genes", "73.88", "73.26", "−0.62 [−1.19, −0.03]"),
          ("DEV, ≥ 6 segments", "71.69", "70.95", "<b>−0.74 [−1.71, +0.23]</b>"),
          ("Set B, pooled", "68.47", "65.06", "−3.41 [−3.78, −3.05]"),
          ("Set B predictions / references", "—", "24,590 / 21,875", "+12% over-calling"),
      ], [36, 16, 22, 26], num_cols=(1, 2))
      + note_box("Criterion, registered before results: d(≥ 6 segments) ≥ +2.58 (twice the T2 seed SD). Verdict: no detectable effect.", RUST)
      + foot("Protocol 33 §8–9 · DEV and set B · 1-Mb block paired bootstrap · one run"),
      """Why M41 was run (protocol 33, registered 11 November). The panel had shown that T2 falls behind on long, many-exon genes. If the grammar's single intron state could not represent both short and very long introns, splitting it into a short-intron and a long-intron state (with their own length behaviour) might help many-exon genes. This was the one grammar-side change we tried for C4.

Criterion (verbatim in protocol 33): the main endpoint is d(≥ 6), the difference in DEV exact-chain F1 on genes with six or more CDS segments between M41 and T2, with the threshold at twice T2's seed standard deviation (2.58 points).

Results (21–22 November, protocol 33 §9). DEV 73.26 overall, 70.95 on ≥ 6-segment genes. Against T2: d(≥ 6) −0.74 [−1.71, +0.23]; overall −0.62 [−1.19, −0.03]. Set B: 65.06, 3.41 below T2 [−3.78, −3.05], with 24,590 predicted chains against 21,875 references (12% over-calling). Verdict: no detectable effect.

What it tells us. The multi-exon and long-gene gap is not fixed by a richer intron length model in the grammar. Together with the long-gene stratification (previous slide), this points to the potentials (backbone features and training geometry) rather than the grammar's capacity, which is what the older C4 replay suggested too. The cost was about one training run on ont; a 17-hour delay came from the TE-gate bug.

Not tried. Training with longer windows so that long genes are seen complete; a decoder that uses evidence across windows. Either would be a new experiment (GPU) and is listed on the questions slide.""")

# ===================================================================================================== dropped
write("dropped", head("Tested and dropped", "Nine directions failed bars we set in advance; each one narrowed the story rather than ending it")
      + table(["Direction", "Test", "Result"], [
          ("Beat Carbon-A with a residual head (M29)", "DEV vs the Carbon-A gate", "71.84 vs 74.85: dropped 17 Oct"),
          ("Adapt a head on frozen Carbon-A (M30)", "U − filtered Carbon-A ≥ +3.0", "−2.32 / −2.73: not met"),
          ("More training species (M34 P2, P2b)", "T10 − T2 ≥ +2.0", "−0.96 DEV, −0.14 set B"),
          ("Self-training on the target genome (M36)", "close ≥ 30% of the gap to Carbon-A", "all 7 arms 1.5–7.7 below zero-shot"),
          ("Joint multi-assembly decoding (M37)", "headroom ≥ 20%", "10.3%: NO-GO"),
          ("Zero-shot near Carbon-A (C1 panel)", "gap ≤ 6.0, ≥ Tiberius, Vipsania", "−12.60: fails"),
          ("Intron-length states (M41)", "≥ 6-segment gain ≥ +2.58", "−0.74: no effect"),
          ("Partial labels on real data (M42)", "U > Z and U > F", "−2.01 and −0.37: stop"),
          ("Add model genes to old annotation", "F1 not lower (post hoc)", "−0.7 to −6.9"),
      ], [40, 30, 30], size=24)
      + foot("Reports 11, 12, 14, 20, 24, 27, 31, 33, 35 · each criterion committed before its result, except the post-hoc merge rule (committed before running)"),
      """Each row: the question, its pre-registered bar, and the number.

M29 (16–17 October): a residual head on Carbon-A's frozen offsets: DEV 71.84 below the Carbon-A gate 74.85 and filtered Carbon-A 77.41; temperature scaling (A2) also failed. We stopped competing on accuracy.
M30 (17–19 October): adaptation on frozen Carbon-A features with 25% labels: U − filtered Carbon-A −2.32 [−2.72, −1.94] / −2.73 (needed ≥ +3.0). The mechanism part (U − F +16.50 / +14.15) held and became C2.
M34 P2 / P2b (26 October): ten species instead of two: DEV −0.96 [−1.52, −0.29], set B −0.14 [−0.53, +0.28]. Arms T4 and T6 were withdrawn; result in Supp. S9.
M36 (31 October): target-genome bootstrap with self-training seeds, protein seeds, their union and a fill control: every arm 1.5–7.7 below zero-shot; closure −26% / −63% (needed +30%); C2 withdrawn for pseudo-labels (F − S1 −0.94 [−1.61, −0.28], under the 2-point bar although F is significantly worse).
M37 (22 October): joint decoding across assemblies: headroom 10.3% < 20%.
C1 panel (17 November): −12.60 against filtered Carbon-A, also below Tiberius (−4.20) and Vipsania ft (−2.36).
M41 (22 November): −0.74 [−1.71, +0.23] against the +2.58 bar.
M42 gate 1 (24 November): U − Z −2.01 [−2.56, −1.50], U − F −0.37 [−0.75, +0.01]; missing-locus recall 15.70 vs 16.15.
Merge (24 November, post hoc with rules committed first): F1 −0.7 to −6.9 against the old annotation.

Earlier (June–September, M1–M27): fixed-caller cross-species transfer, CRF and Tversky losses on per-base models, the M25R structure head and many M28 candidate pools also failed their bars (report 03).

Why keep them. Each is pre-registered, so they are honest negative results. Two of them (M36, M42) are directly relevant to C2's limits and belong in the paper; most others belong in the supplement.""")

# ===================================================================================================== integrity
patch("integ-stats", [
    ("Every number is one training run; the CIs resample genes or 1-Mb blocks, never seeds", "Most numbers are one run; T2 and P4 now have seeds, and their spread (up to 3 points) matches several effects"),
    ("One training run per model; two more T2 seeds queued", "T2: three runs (DEV 71.65–74.09); P4: two per arm; the rest one run"),
    ("dropped for our models, not for the tools", "dropped for every model and tool on the panel"),
])
patch("integ-exposure", [("and only the panel is fresh; even it contains set B’s two species", "and the panel was fresh when scored; it contains set B’s two species")])

# ===================================================================================================== synthesis
write("synthesis", head("Synthesis", "What the evidence supports: a safe way to train on incomplete labels and a useful gene confidence score, not a more accurate annotator")
      + '<div style="display:grid;grid-template-columns:1fr 1fr 1fr;gap:20px">\n'
      + "".join(f'<div style="display:flex;flex-direction:column;gap:12px;background:{c};border-radius:8px;padding:22px 26px">\n<h3 style="font-size:30px;font-weight:600;color:{INK}">{t}</h3>\n'
                + "".join(f'<p style="font-size:24px;line-height:1.35;color:{BODY}">{x}</p>\n' for x in xs) + "</div>\n"
                for t, c, xs in [
                    ("Holds", "#E4ECF3", ["One grammar layer trains on three backbones; the backbone sets accuracy (+8.6 PlantCAD2 over NT-v2)",
                                          "Posteriors rank genes (AUROC 0.77); a frozen threshold helps 10 of 10 species"]),
                    ("Holds with limits", PANEL, ["Marginalising beats filling when missing genes are known (M30; P4 in 2 seeds), but training is unstable",
                                                  "Posteriors are over-confident across species"]),
                    ("Failed", "#F5E6DC", ["Zero-shot near Carbon-A (−12.6 on the panel)", "Partial labels on a real incomplete annotation (M42)",
                                           "More species, self-training, intron states"]),
                ])
      + "</div>\n"
      + note_box("<b style=\"color:#1B2430\">The shared weakness:</b> long genes. Above 10 kb T2 recalls 0–24% of genes against Carbon-A’s 46–59%, mostly by splitting them.")
      + foot("All claims slides · reports 20, 21, 31, 35, 36, 37"),
      """How the pieces fit. The grammar layer gives two properties that the competing tools do not have by construction: training with unknown positions summed out, and gene-level posteriors from the partition function. Accuracy, on the other hand, comes from the backbone: swapping NT-v2 for PlantCAD2 under the same layer moved DEV by +8.6 and set B by about +20 points, while adding species, self-training or a richer intron model did not help.

The properties hold in a narrow form. Marginalising is safer than filling whenever we know which genes are unlabelled (M30, P4 two seeds); but the oracle-free test M42 shows that, with a real old annotation and protein evidence, training on the target genome does not beat the zero-shot model, and M36 showed the same with pseudo-labels. The posterior is a good ranking score and a portable filter, but not a calibrated probability across species.

The accuracy gap to Carbon-A on the independent panel (12.6 points median) is concentrated in long genes, which T2 splits. That is the clearest target for future work and the most likely reviewer question.

What this means for the story (to discuss on the questions slide). The honest paper is about how to learn gene structure from the incomplete annotations that real genomes have, and what that buys and does not buy: a method plus an honest benchmark. It is not a new state of the art.""")

# ===================================================================================================== status
write("status", head("Status · 24 November 2026", "All planned analyses are done; insects are paused for TE until early December, and GeneCAD waits for your budget decision")
      + cards([
          ("Done", "Panel and C1 verdict", "C1 · Fig. 4 · scored 17 Nov", "Six tools on ten plants; C1 fails (−12.60); seed rule triggered, not run"),
          ("Done", "M42 gates 0–1, merge", "C2 · 24 Nov · 41.5 yuan", "Gate 1 stops (U − Z −2.01); merge lowers F1; route A"),
          ("Done", "C3 and C4 final", "Reports 36, 37 · 24 Nov", "AUROC 0.77, over-confident; long genes split; P4 unstable"),
          ("Paused", "M39 insect replication", "C2 · NT-v2 · protocol 29", "I10_k50U at step 2,304; TE has ont since 23 Nov; result about 3 Dec"),
          ("Needs your decision", "GeneCAD baseline", "C1 · Table 1 · closest method", "CAU about 240–340 yuan, or Baobab after 1 Dec"),
          ("Not started", "Panel seeds, M40, writing", "Spread only · 0% endpoint", "Panel seeds after 1 Dec; M40 last; writing waits for the story"),
      ])
      + foot("Reports 20, 29, 31, 35, 36, 37 · dashboard of 24 November"),
      """Card by card.

1. Panel and C1 (done 17 November). Core-10 panel, 119 sequences, 4.79 Gb; T2 inference 33.8 GPU hours (3 November); six baselines about 85 GPU hours on ont (9–15 November). C1 fails. The seed rule registered on 7 November was triggered (T2_s2, T2_s3 on the panel) but has not run; it only adds spread.

2. M42 (done 24 November). Gate 0 on the Mac, gate 1 on CAU (U, F, U3, Z; about 16.85 GPU hours, 41.5 yuan with smoke and builds). Stop. The merge analysis lowers F1. Your decision: route A plus the merge analysis.

3. C3 and C4 (done 24 November). Reports 36 and 37, plans committed before computing, all on the Mac from copies of ont's outputs (no ont CPU used).

4. M39 (paused). The insect replication of P4 on NT-v2 (ten insects, P4-I: k50U − k50F ≥ +2.0 pooled and ≥ 6 of 8 species). I10_k50U reached step 2,304 (selection 58.69 / 34.02 / 59.28 at 768 / 1,536 / 2,304). TE's timing benchmark has the GPU since 23 November; it will tell us when it resumes (about 29 November to 2 December); I10_k50F and I10_T follow, result about 3 December.

5. GeneCAD (your decision). Installable on CAU (Python 3.12, torch 2.7.1, flash-attn and mamba compiled from source). Published speed about 18 kb/s on an A100; on a 4090 an estimated 10–14 kb/s, so 100–140 GPU hours for the 5.1 Gb panel, about 240–340 yuan. Alternative: Baobab after 1 December (asked first).

6. Not started. Panel seeds (after 1 December). M40, the semi-supervised 0%-label endpoint (code and CPU tests only). Writing: the story architecture of 24 November (commit 5114242) is a proposal; drafting waits until we settle the story together.""", bg=PANEL)

# ===================================================================================================== branches
write("branches", head("Decision branches · pending results", "Two pending results can still change the paper; both have their readings written in advance")
      + table(["Pending", "If it holds", "If it fails"], [
          ("<b>M39 insects</b> (P4-I: U − F ≥ +2.0, ≥ 6 of 8)", "“The C2 advantage replicates in insects (NT-v2, one seed)”: C2 spans two lineages and two backbones", "C2 stays plant-only; cannot tell backbone from lineage without an NT-v2 plant pair (≈ 50 GPU h)"),
          ("<b>GeneCAD</b> on the panel (no criterion)", "If close to Carbon-A: our gap is about the layer, not PlantCAD2", "If close to T2: the gap is the backbone; table complete either way"),
          ("<b>Panel seeds</b> (spread only)", "Report the spread", "C1 verdict unchanged"),
          ("<b>The story</b> (your decision)", "Method + honest benchmark (Genome Biology)", "Analysis paper, or a narrower method paper (NAR GB)"),
      ], [26, 37, 37], size=24)
      + foot("Protocol 29 §4.1 (P4-I, registered before the insect runs) · report 28 §7 · GeneCAD: no criterion registered"),
      """M39 (protocol 29 §4.1, registered before the insect runs). P4-I: on the eight new insects' validation regions, k50U − k50F ≥ +2.0 pooled and the same sign in ≥ 6 of 8 (≥ 3/4 if species had to be replaced). Wording under each outcome, as registered: both conditions met → "the C2 advantage replicates in insects (NT-v2, one seed)"; 0 < pooled < +2.0 with enough species → "same direction, below the bar"; otherwise not replicated. Note in the protocol: insects change both the backbone (NT-v2) and the lineage, so a failure cannot be attributed without an NT-v2 plant pair (about 50 GPU hours, not planned). Also: the M39 insect arms use the same oracle mask as P4, so even a success stays an upper-bound result.

GeneCAD. No criterion has been written, and it should be before the result: I suggest it is descriptive only (a row in the panel table and the long-gene stratification), because it shares our backbone (PlantCAD2) and adds a different head and a chromosome-wide CRF. Readings I would write: if GeneCAD is near Carbon-A, our gap is in the grammar layer or training, not the backbone; if near T2, the backbone limits both.

Panel seeds (report 28 §7, registered 7 November): spread only; the C1 verdict stays as scored.

The story is yours to decide (next slide). The writing skill's rule: no core claim, no article. I think we still have one, in narrower form; the question is how to frame it for Genome Biology.""")

# ===================================================================================================== timeline
write("timeline", head("Timeline", "No submission date yet: the next fixed points are ont returning from TE (about 1 December) and the Baobab reset on 1 December")
      + table(["When", "What", "Depends on"], [
          ("This week", "Review this deck together; decide the story, the claims and what (if anything) to add", "you"),
          ("Now", "GeneCAD on CAU (if approved) · figure and Results drafting once the story is set", "your budget decision"),
          ("≈ 29 Nov–2 Dec", "ont back from TE: M39 resumes (I10_k50U, then I10_k50F, I10_T)", "TE timing benchmark"),
          ("1 Dec", "Baobab credits reset: panel seeds, GeneCAD if not on CAU (asked first)", "credit reset"),
          ("≈ 3 Dec", "M39 verdict; C2 one or two lineages", "M39"),
          ("Then", "Methods, Discussion, Introduction (facts checked with ChatGPT), abstract", "story, M39"),
          ("Then", "Mock review: ChatGPT, AI Studio, Claude on the main text only", "full draft"),
      ], [16, 60, 24], size=24)
      + foot("Writing skill (story-driven-manuscript) order · dashboard of 24 November · no date is promised"),
      """The writing workflow you set (story-driven-manuscript): results complete and figures done, then agree the story and the journal, then write in human order: figures and Results, Methods, Discussion, Introduction, abstract and title last. Every round gives an English PDF, a Chinese review version, a list of changes and the decisions I need from you, then stops for your feedback.

Fixed dates. ont returns from TE's timing benchmark when the TE session says so (it will send one line); its estimate was 29 November to 2 December. Baobab's billing year resets on 1 December; each job there is asked first.

What could still move. GeneCAD (if you approve CAU, it could run this week: 100–140 GPU hours on up to two cards, about 3 days). M39 (about 3 December). A long-gene experiment, if you want one (training with longer windows): a new GPU experiment, not scheduled.

No submission date is promised: it depends on the story decision and on whether a new experiment is added.""")

# ===================================================================================================== limits
write("limits", head("Limits", "Six limits any reviewer will raise, each with its number")
      + cards([
          ("Few runs", "P4 margin +3.24 vs +14.97 between seeds; T2 runs span 3 points; most models one run"),
          ("Upper-bound setting", "C2’s tests placed the unknown mask with the full truth; the real test (M42) failed"),
          ("Behind the best tools", "Panel median 12.6 below filtered Carbon-A, 4.2 below Tiberius"),
          ("Long genes", "Recall 24%, 3% and 0% for 10–20, 20–50 and ≥ 50 kb, against 59, 54 and 46%"),
          ("Over-confident posteriors", "P ≥ 0.9 bin right 77% pooled, 58% in Themeda"),
          ("Computed truth, small compute", "RefSeq or submitter models, no second truth; one shared GPU"),
      ], badge=False, size=26)
      + foot("Reports 20, 31, 35, 36, 37 · integrity slides"),
      """1. Runs. T2 three runs: DEV 73.88 / 74.09 / 71.65, new-species regions 68.84 / 67.28 / 65.83, set B 68.47 / 66.86 / 65.93. P4 two runs per arm with very different margins. Everything else one run. Seed spread (1.3–1.5 SD) is as large as several effects (T10 − T2 −0.96, M41 −0.74).
2. Upper bound. M30 and P4 marked exactly the removed genes as unknown. M42 replaced that with protein evidence on a real old annotation and stopped; M36 with pseudo-labels failed too.
3. Accuracy. Panel medians: T2 55.05, Tiberius 59.25, Vipsania ft 57.41, ANNEVO 63.08, filtered Carbon-A 67.65. On submitter-annotated species T2 is level with Tiberius and Vipsania (56.6 vs 53.3 and 56.4).
4. Long genes. 8.4% of panel genes, half of the recall gap; splitting dominates in Hakea, Amborella, Akebia.
5. Calibration. Pooled ECE 0.088; five species miss the high-confidence bin by more than 0.15.
6. Truth and compute. The truth is one transcript per gene from RefSeq (Gnomon) or submitters, no RNA-seq or second truth; DEV is in PlantCAD2 pretraining and Carbon-A training. One shared RTX 4090; panel baselines run once.""")

# ===================================================================================================== questions
write("questions", head("Questions to settle together", "Before writing, we need to agree on the story; these are the choices I see, with my view")
      + table(["Question", "Options", "My view"], [
          ("1. What is the paper about?", "(a) a method: training on incomplete labels + gene confidence, with an honest benchmark; (b) an analysis: what limits DNA-LM gene annotators (backbone, labels, length)", "(a), with (b)’s findings as its second half"),
          ("2. Is Genome Biology right?", "GB · NAR GB · Bioinformatics", "GB is a stretch without a biological result; NAR GB fits"),
          ("3. Add an experiment?", "GeneCAD (240–340 yuan) · long-window training for long genes · nothing", "GeneCAD yes; long windows only if it is the headline"),
          ("4. Where do negatives go?", "M42 and M36 in main text (C2’s limits); M41, M34, M29 in supplement", "as listed"),
      ], [22, 50, 28], size=24)
      + foot("Writing skill rules: no core claim, no article; journal, title and authors are your decisions"),
      """1. The story. Option (a), a method paper: "a probabilistic gene grammar lets DNA language models learn from incomplete annotations without being taught that missing genes are background, and gives a gene-level confidence; on independent plants it is honest about where it falls short (long genes)". Option (b), an analysis paper: "what limits DNA language model gene annotators: the backbone sets accuracy, incomplete labels hurt unless marginalised but do not add information, and long genes break". My view: (a) with the findings of (b) as its second half; it matches the evidence and keeps C2 as the core claim in its narrow form.

2. Journal. You chose Genome Biology on 24 November. Genome Biology editors usually want a biological result or a broadly useful tool; we have neither a state-of-the-art tool nor a new biological finding. NAR Genomics and Bioinformatics fits a method-plus-benchmark paper with negative results. I would keep GB only if we add something that changes the story (question 3).

3. Experiments. GeneCAD makes the benchmark table complete and answers the most likely reviewer question (same backbone, different head); about 240–340 yuan on CAU or free on Baobab after 1 December. A long-window training run could test the long-gene hypothesis and, if it worked, would be a real result; it is new GPU work (several days) and would need its own protocol. M39 runs anyway when ont returns.

4. Negatives. M42 and M36 define C2's limits and belong in the main text; M41, the multi-species result, M29/M30 accuracy attempts and M37 go to the supplement as pre-registered negatives.

After we settle these, I will write the story architecture again and start with figures and Results.""")

# ===================================================================================================== deck.json
removed = ["c1-scale", "c2-ckpt", "c3-pc2", "branches2"]
for r in removed:
    p = os.path.join(S, f"{r}.html")
    if os.path.exists(p):
        os.remove(p)
order = ["cover", "problem", "journey", "problems", "changes", "claims", "method",
         "c1-dev", "c1-setb", "c1-panelset", "c1-panel",
         "c2-m30", "c2-p4", "c2-dyn", "c2-m36", "c2-critique", "c2-m42design", "c2-m42", "c2-merge",
         "c3-calib", "c3-rel", "c3-sel",
         "c4-strata", "c4-replay", "c4-long", "c4-m41",
         "dropped", "integ-exposure", "integ-truth", "integ-stats",
         "synthesis", "status", "branches", "timeline", "limits", "questions"]
for sid in order:
    assert os.path.exists(os.path.join(S, f"{sid}.html")), sid
d = json.load(open(os.path.join(P, "deck.json")))
d["title"] = "coding-rna 系统梳理 · 2026-11-24"
d["order"] = order
d["sections"] = {
    "s0": {"description": "The problem, how we got here, the problems we hit, what changed, the claims and the model", "start": "cover"},
    "s1": {"description": "C1: one grammar layer on several backbones; the panel verdict", "start": "c1-dev"},
    "s2": {"description": "C2: training on incomplete annotation, its critiques and the oracle-free test", "start": "c2-m30"},
    "s3": {"description": "C3: gene posteriors on the final model", "start": "c3-calib"},
    "s4": {"description": "C4: where the errors are; long genes", "start": "c4-strata"},
    "s5": {"description": "What was tested and dropped, and how far the evidence can be trusted", "start": "dropped"},
    "s6": {"description": "Synthesis, status, pending results, timeline, limits and the questions to settle", "start": "synthesis"},
}
json.dump(d, open(os.path.join(P, "deck.json"), "w"), indent=2, ensure_ascii=False)
print("ok", len(order), "slides")
