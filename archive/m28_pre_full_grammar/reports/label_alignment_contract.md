# M28 R4 target alignment — before fitting

Use the frozen M26/M28 CDS-assessable policy for choosing a complete primary:
longest assessable CDS per gene, tie transcript ID; then apply training-only
support/overlap masks. Unsupported but assessable primaries stay in DEV truth.
For a gene with no assessable isoform, retain the old longest-all local-only
primary to preserve certain grouped CDS/splice supervision; it is never a full
chain positive. Other isoforms and unknown structures retain negative protection.

Two CPU jobs, each 2CPU/8GiB/10min, train/DEV allowlist only, zero optimizer steps.
New outputs only. Reconstruct targets/unknown intervals on the unchanged natural
grid; assert DEV assessable primary IDs exactly match the frozen ruler.

Do not resample: keep the actual R2 paired 1536 draws, their order and original
q/importance weights. The earlier 75% uniform component gives positive support
to every natural-grid window; weighting still targets uniform windows within
each fixed stratum after label changes. Do not recompute q under a law that did
not generate these draws. Record changed targets and actual unique exposure.
Feature tensors do not change; trainers must join updated labels by window ID,
not take obsolete supervision metadata from R3 feature cache objects.
