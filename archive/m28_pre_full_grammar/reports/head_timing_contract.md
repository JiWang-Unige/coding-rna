# M28 R3 head timing — fixed before execution

Use the same twelve preselected TRAIN feature caches from the successful fixed
GENERanno revision smoke. No GLM re-extraction, validation/test/Setaria reads,
reference-dependent example selection, optimizer or parameter update.

Each arm C0/B1 uses seed 0, FP32 trainable heads, hidden 128, identical feature,
DNA, window order, labels and per-window sampling weights. One forward and
backward per cached window/arm. C0 uses the shared local fine/grouped objective.
B1 adds masked class-balanced endpoint/link/chain losses (all coefficients 1).
Generate and save free B1 candidates before reference labels or injection.
Candidate budget is unchanged from integration. This is throughput and actual
gradient plumbing, not a head fit or predictive performance comparison.

Record feature I/O, shared forward, free proposal, supervised head construction
plus backward, parameter counts and peak memory separately. Confirm actual
positive, reference-negative and unknown link-logit gradient directions; require
positive and negative connection supervision in at least one existing cache for
each species. No candidate/sample selection by performance. No throughput-based
accuracy claim; isolated windows are not a whole-chromosome cost guarantee.

One shared RTX3090/2CPU/16GiB job, 15min ceiling (0.25 GPUh); no automatic retry or
expanded sweep. Freeze a separate bounded training contract after the measurement.
