"""Record the CRITICAL SegmentNT segmentation-leakage finding (M8-CK4 audit, surfaced early).
Appends docs/10 finding + annotates docs/06 M7 result + docs/00. python3 scripts/_log_segnt_finding.py"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

D10 = r"""

## Finding (Research/Comparability) 2026-06-12 -- CRITICAL: SegmentNT(multi_species) was SEGMENTATION-FINE-TUNED on chicken+fly -> our chicken/fly evals are LABEL-LEAKAGE contaminated; only plants (arabidopsis) are truly clean held-out (M8-CK4 audit, surfaced early)
SegmentNT(segment_nt_multi_species) audit (HF model card + dataset card, WebFetch-verified): backbone = nucleotide-transformer-v2-500m-multi-species (850 genomes, 174B nt) which EXPLICITLY EXCLUDES plants+viruses; segmentation head FINE-TUNED on human + 5 species = {mouse, chicken, fly, zebrafish, worm}. CONSEQUENCE for our project (candidate = frozen SegmentNT 14-element features as input):
- Gallus gallus (chicken): SegmentNT segmentation head was SUPERVISED-TRAINED on chicken gene annotation -> feeding chicken SegmentNT features ~= reading off a model that already learned chicken gene structure. The M7 chicken result (candidate spec 0.970 vs raw-DNA anchor 0.669, +0.30) is LABEL-LEAKAGE inflated, NOT clean cross-clade generalization.
- Drosophila melanogaster (fly): ALSO in the 5 segmentation species -> the ORIGINAL FP-FRAGFIX-CONSTR yeast+fly result (spec 0.9218) had the fly half contaminated too.
- Arabidopsis thaliana: NEITHER backbone (plants excluded) NOR segmentation (no plant) saw it -> TRULY CLEAN held-out. The candidate STILL beats the anchor on arabidopsis (spec 0.954 vs 0.892, +0.06) -> this is the clean, real generalization signal (and on a kingdom SegmentNT never saw).
- Saccharomyces (yeast, fungi): segmentation-clean (not in the 5); backbone may include fungi bucket -> weaker (representational) overlap only.
LESSON: the held-out re-anchor's headline (+0.155 spec cross-clade) is INFLATED by chicken leakage; the CLEAN cross-clade margin is the arabidopsis +0.06. SegmentNT-feature gains are substantially confounded by "SegmentNT already did gene-segmentation on this exact species" for chicken/fly. PRE-CLAIM GATE: any full/scale claim must evaluate on species NOT in {human,mouse,chicken,fly,zebrafish,worm} (segmentation-clean) and ideally not in the backbone (plants). M8 redirected: clean held-out = {arabidopsis, rice (Oryza sativa)} (both-layer clean plants); chicken micro+macro retained as a CONTAMINATED robustness stratum (labeled, not held-out generalization). This is exactly the leakage the retrospective + tri-review pre-claim guard warned about — now confirmed decisively at the segmentation-label level (worse than the suspected backbone-only overlap).
"""

D06_ANNOT = """

### CRITICAL POST-HOC ANNOTATION (2026-06-12, M8-CK4 SegmentNT audit) — REANCHOR-HELDOUT-M7 chicken result is LEAKAGE-contaminated
SegmentNT(multi_species) segmentation head was FINE-TUNED on chicken (one of {human,mouse,chicken,fly,zebrafish,worm}). So the M7 chicken candidate result (per-species spec 0.970 vs anchor 0.669) is LABEL-LEAKAGE inflated — SegmentNT already learned chicken gene structure. ONLY arabidopsis (plant; excluded from BOTH backbone and segmentation training) is a TRULY CLEAN held-out: there the candidate still wins (spec 0.954 vs 0.892, +0.06) — that is the honest cross-clade signal. The base-weighted M7 headline spec 0.9604 (+0.155 over anchor) is INFLATED by the contaminated chicken half; clean margin ≈ +0.06 (arabidopsis only). Fly (in the original yeast+fly eval) is also contaminated. See docs/10 (2026-06-12 SegmentNT finding). M8 redirected to add clean species (rice) + flag chicken as contaminated-robustness-only.
"""

def app(p, t):
    open(os.path.join(ROOT, p), "a").write(t); print("appended", p)

app("docs/10_findings.md", D10)
app("docs/06_results_log.md", D06_ANNOT)

p = os.path.join(ROOT, "docs/00_active_goal.md"); t = open(p).read()
mk = "## last_result_summary\n"; i = t.index(mk) + len(mk)
blk = (
"- **CRITICAL 2026-06-12 (M8-CK4 SegmentNT audit, surfaced early)**: SegmentNT(multi_species) segmentation head was FINE-TUNED on {human,mouse,chicken,fly,zebrafish,worm} -> our **chicken + fly** evals are LABEL-LEAKAGE contaminated; **only arabidopsis (plant) is truly-clean held-out** (candidate still wins there: spec 0.954 vs 0.892, +0.06 = the honest signal). M7 headline +0.155 is inflated by chicken leakage. PRE-CLAIM GATE: claim only on segmentation-clean species (not the 6). M8 redirected: clean held-out = {arabidopsis, rice}; chicken = contaminated robustness stratum (labeled). See docs/10.\n")
open(p, "w").write(t[:i] + blk + t[i:]); print("docs/00 annotated")
