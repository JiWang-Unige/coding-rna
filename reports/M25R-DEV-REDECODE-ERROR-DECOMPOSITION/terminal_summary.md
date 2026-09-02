# M25R Stage 1 terminal summary

- Analysis: `M25R-DEV-REDECODE-ERROR-DECOMPOSITION`
- Slurm job: `12128521`
- Terminal state: `FAILED`, exit `1:0`, elapsed `1-00:12:36`
- Experiment status: `FAILED`
- Code commit used by the run: `1641ff72e3e68fb4dfb91f1fc228322113e154c7`

Two earlier attempts are preserved on Baobab but are not scientifically interpretable: job `12127330` exited after `14:34` at `384/1,536` training windows, and job `12127541` exited after `36:08` at `1,152/1,536`. Their tracebacks were truncated before the exception body. R2 changed only Python log placement to node-local storage and produced the complete terminal traceback below; it does not retrospectively prove the precise causes of those earlier exits.

## Frozen scope

The run used only the three existing M25R checkpoints and the frozen highest-ranked decoder tuple for each epoch. It used all `1,536` observed M25R training windows and the two Arabidopsis/rice development validation chromosomes. It did not update weights, search thresholds or decoder settings, or read Setaria files. `resolved_inputs.json` records `setaria_files_read=false`, `threshold_or_decoder_search=false`, and `weights_updated=false`.

## What completed before termination

Epoch 1 completed raw-head inference on all `1,536` training windows and inference on both development validation chromosomes. Control flow then passed:

1. reproduction of the four frozen aggregate metrics within absolute tolerance `1e-5`;
2. assignment of all `6,450` reference chains to one earliest attrition stage;
3. writing of the epoch-1 attrition table, candidate-lineage ledger and replayed GFF3;
4. reconciliation of `2,098` emitted models with the GFF3 gene/transcript counts and `100%` audit coverage.

The exact replayed aggregate values were not persisted before termination, so only passage of the assertions—not the values themselves—is established by the control flow.

## Terminal failure

The independent GFF3 audit checked all `2,098` emitted transcripts. Only `1,338` were structurally valid; `760` failed, giving validity `0.6377502383`. The script then raised at `redecode_error_decomposition.py:988`:

```text
AssertionError: independent GFF3 validity audit failed
```

This is not a training crash and not a Slurm, GPU-memory or logging failure. The batch step used about `9.1 GiB` maximum RSS out of `96 GiB`. The earlier node-local logging change succeeded and preserved the complete traceback.

The immediate scientific result is that the frozen epoch-1 decoder emits many structurally invalid models. The immediate implementation defect is that the diagnostic treats the phenomenon it was designed to measure as a fatal integrity error. It therefore stops before serializing the epoch result, before running epochs 2 and 3, and before writing `stage1_diagnostic.json`.

## Information still missing

- No final three-epoch diagnostic exists.
- Epochs 2 and 3 were not evaluated.
- The `760` invalid transcripts are labeled only `structural_validity`; the saved exception does not separate phase continuity, start/stop, frame length, internal stop, splice motif or codon-feature mismatch.
- `audit_gff3()` returns only the first 50 generic failures and does not persist a complete component-level validity ledger.
- Stage-level attrition and error summaries were computed for epoch 1 but were not written before the fatal assertion.

These are missing measurements, not zero-valued results.

## Code-review question

Review the following files together:

- `scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/redecode_error_decomposition.py`, especially `audit_gff3()` and `run_epoch()`;
- `scripts/eval_m25_structure.py`, especially `structurally_valid()`;
- `tests/test_m25r_redecode_error_decomposition.py`;
- `sbatch/M25R-DEV-REDECODE-ERROR-DECOMPOSITION.sbatch`;
- this directory's `python_12128521.out`, `python_12128521.err`, and `resolved_inputs.json`.

The review must decide the smallest way to make diagnostic failures durable without weakening integrity checks. In particular:

1. Should a full-coverage validity audit with invalid predictions be recorded as a scientific result while audit incompleteness or count mismatch remains fatal?
2. Should each epoch write an atomic result before the next epoch begins so one failed epoch cannot erase completed evidence?
3. Which component-level validity outcomes must be recorded for every transcript?
4. How should the script continue through all three frozen epochs while still ending in a review-required, non-success scientific status?
5. Which focused tests would have caught the current fail-fast behavior before a 24-hour cluster run?

No code repair or rerun is authorized by this summary. The next action is review.

## Archived compact artifacts

- `JOBID`
- `STATUS`
- `env.txt`
- `resolved_inputs.json`
- `python_12128521.out`
- `python_12128521.err`

The 14 MB candidate-lineage ledger, replayed GFF3, attrition TSV, model checkpoints and all source data remain on Baobab and are intentionally excluded from Git.
