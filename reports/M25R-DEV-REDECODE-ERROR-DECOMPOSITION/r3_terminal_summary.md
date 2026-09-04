# M25R Stage 1 R3 terminal failure summary

- Analysis: `M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R3`
- Slurm job: `12166448`
- Code commit: `05b4f4e31ba6b055cf46c39b89e28002032d4142`
- Terminal state: `FAILED`, exit `1:0`, elapsed `1-01:03:39`
- Peak RSS: `9,570,072 KiB`
- Experiment status: `FAILED`

## Frozen scope

The job used the three existing M25R checkpoints and the frozen highest-ranked epoch tuple on the Arabidopsis/rice development primary nuclear chromosomes. It did not update weights, search thresholds, alter the decoder or read Setaria files. `resolved_inputs.json` records `setaria_files_read=false`, `weights_updated=false` and `threshold_or_decoder_search=false`.

## Evidence completed and preserved

Epoch 1 completed raw-head inference on all 1,536 frozen training windows and validation inference on `arabidopsis_thaliana/NC_003074.8` and `oryza_sativa/NC_089041.1`. Reaching the later failure proves that the frozen aggregate reproduction assertion passed within `1e-5`.

The following epoch-1 artifacts were persisted before termination:

- `6,450`-row `reference_attrition.tsv`;
- `27,913`-row `candidate_lineages.jsonl`;
- `2,098`-model `replayed_predictions.gff3`;
- `2,098`-row `structural_validity.jsonl`.

Reference-chain earliest-stage counts were: region state path `105`, no covering non-intergenic block `1,922`, wrong CDS-run count `586`, start/stop motif candidates `601`, ordered donor/acceptor candidates `473`, learned boundary choice `229`, phase check `1,145`, complete-ORF/internal-stop check `1`, boundary threshold filter `2`, and emitted exact chain `1,386`. No reference was assigned to legal-terminal-transition failure.

Fixed post-hoc counts over all 6,450 references were: transition reachable `3,727` (`0.57783`), motif reachable `1,693` (`0.26248`), and truth-assisted exact chain `895` (`0.13876`). These are diagnostics, not model-performance claims.

The independent GFF3 ledger reproduced the prior epoch-1 finding: `1,338/2,098` transcripts valid (`0.6377502383`) and `760` invalid. Component failures were start-codon sequence `417`, stop-codon sequence `309`, and splice motif `167`; counts overlap across transcripts. Parent linkage, phase values/continuity, codon feature rows, minimum length, frame length and internal-stop checks had zero failures.

## Terminal failure and root cause

After writing the complete structural-validity ledger, Python raised:

```text
KeyError: 'epoch'
```

`run_epoch()` receives the frozen validation-grid tuple in parameter `row`. The new ledger writer used `for row in transcript_ledger`, and Python function scope therefore replaced the frozen tuple with the last ledger record. Construction of `frozen_tuple` later evaluated `row['epoch']` and failed.

This is a deterministic local-variable shadowing defect, not a Slurm, GPU, memory or scientific-model failure. Epoch 1 `diagnostic.json` was not written, epochs 2 and 3 were not started, and no final Stage 1 report exists. The prior unit test monkeypatched `run_epoch()` and the CPU smoke exercised `audit_gff3()` only, so neither test covered post-audit result assembly.

## Interpretation boundary

The persisted epoch-1 ledgers are valid partial evidence, but Stage 1 is incomplete and cannot yet support a final representation/supervision/decoder/data-transfer go/no-go. The result continues to show severe upstream attrition and poor structural validity, while Setaria remains unaccessed.

## Review-required next decision

The minimal proposed repair is to rename the ledger loop variable, add one focused test that executes post-audit result assembly with a real frozen-tuple-shaped row, and use a new immutable output ID if a rerun is approved. Review must decide whether that is sufficient and what exact test must pass before another multi-day run.

No code repair or rerun is authorized by this summary.
