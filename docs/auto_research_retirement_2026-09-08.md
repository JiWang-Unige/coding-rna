# Auto Research retirement — 2026-09-08

User requested a conversation-driven research workflow without machine goals or mandatory third-party review.

Recoverable archive: `/home/users/j/jwang/retired-auto-research/coding-rna-20260908T084106Z` (owner-only directory). No permanent deletion.

Removed active project skill mirrors, hooks, startup templates, installation/regeneration scripts, machine goal, MCP configuration and deterministic framework gates. Removed cli_review and architecture-selection policy from cluster_config.yaml; retained compute fields.

Conflicts observed: RUN_PROMPT required absent AGENTS.md; hooks re-injected retired instructions; context_pack claimed exclusive authority; goal signatures caused repeated approval; reviewer retry/quorum rules differed across Skill/config; pre_submit_gate required review records even after conversational approval.

38 historical framework-dependent sbatch launchers were archived unchanged rather than left broken or silently modified. Their source/config/results and scientific protocols remain. The archive is for recovery/provenance, not a runnable standalone checkout.

Preserved: goals/*.md historical protocols, docs/reports/refs/wiki, all data/outputs/checkpoints, research environments, shared account tools/config and other projects. Existing r4_terminal_review.md is untouched. No jobs submitted/cancelled; no Git commit or push.

Local instructions plus the coding-rna-baobab skill are sufficient for local-driven SSH work. No remote AGENTS.md replacement was created. README.md states the new workflow and retires generic historical review/goal requirements without relaxing experiment-specific scientific controls.

## Archived paths

- `.agents`
- `.claude`
- `.codex`
- `agents`
- `.spec-workflow`
- `.mcp.json`
- `.mcp.json.example`
- `ACTIVE_GOAL.json`
- `ACTIVE_GOAL.json.example`
- `CLAUDE.md`
- `CLAUDE.md.example`
- `RUN_PROMPT.codex.md`
- `README.auto-research.md`
- `install.sh`
- `scripts/hooks`
- `scripts/context_pack.py`
- `scripts/iter_ledger.py`
- `scripts/pre_submit_gate.py`
- `scripts/research_flow_guard.py`
- `scripts/repair_advisor.py`
- `scripts/validate_goal.py`
- `scripts/validate_stage_c.py`
- `scripts/build_codex_skills.py`
- `scripts/sync_agents_md.sh`
- `scripts/validate_codex_skills.py`
- `tests/test_validate_goal_profiles.py`
- `README.auto-research.md.backup-20260610-102213`
- `.mcp.json.example.backup-20260610-102213`
- `agents.backup-20260614-214655`
- `scripts.backup-20260610-102213`
- `ARCHITECTURE.md.backup-20260614-214655`
- `CLAUDE.md.framework-upgrade-backup-20260614-214748`
- `.agents.backup-20260614-214655`
- `.codex.backup-20260614-214655`
- `README.auto-research.md.backup-20260614-214655`
- `ACTIVE_GOAL.json.bak-20260614`
- `ARCHITECTURE.md.backup-20260610-102213`
- `.claude.backup-20260614-214655`
- `.mcp.json.example.backup-20260614-214655`
- `RUN_PROMPT.codex.md.backup-20260614-214655`
- `sbatch/M20-STRUCTURED-DECODER-IMPL-SMOKE3.sbatch`
- `sbatch/M22-GENERANNO-1P2B-NONCRF-GBTVERSKY.sbatch`
- `sbatch/M21-GENERANNO-1P2B-CRF-SCREEN-s1-opt.sbatch`
- `sbatch/M18-MULTICLADE-TRAIN-DIAGNOSTIC.sbatch`
- `sbatch/M21-GENERANNO-1P2B-CRF-DIAG-s1.sbatch`
- `sbatch/M16-MULTISPECIES-TRAIN-DISTANCE-SCAN.sbatch`
- `sbatch/M22-GENERANNO-1P2B-NONCRF-GBTVERSKY-SMOKE.sbatch`
- `sbatch/M15-GENERANNO-LORA-PANEL-SCREEN.sbatch`
- `sbatch/M17-SAMEPANEL-GENERALIZATION-BASELINES.sbatch`
- `sbatch/M11-L12-SPEC-CALIBRATION.sbatch`
- `sbatch/M11-L12-SPEC-CALIBRATION-SMOKE.sbatch`
- `sbatch/M20-STRUCTURED-DECODER-IMPL-SMOKE2.sbatch`
- `sbatch/M10-M9L12-CLEANPLANTS.sbatch`
- `sbatch/M12A-FIXEDMODEL-CROSSSPECIES.sbatch`
- `sbatch/M21-GENERANNO-1P2B-CRF-SCREEN-s1-fastval.sbatch`
- `sbatch/M10-GENERANNO-LORA-3C-SMOKE.sbatch`
- `sbatch/M21-GENERANNO-1P2B-CRF-SCREEN.sbatch`
- `sbatch/M20-STRUCTURED-DECODER-IMPL-SMOKE.sbatch`
- `sbatch/M10-GENERANNO-LORA-3C.sbatch`
- `sbatch/M18-GENERANNO-0P5B-SPEC-OBJECTIVE.sbatch`
- `sbatch/M18-GENERANNO-1P2B-SPEC-OBJECTIVE.sbatch`
- `sbatch/M22-GENERANNO-1P2B-NONCRF-GBTVERSKY-s0.sbatch`
- `sbatch/M23-NTV2-CLEAN-TRANSFER-s0.sbatch`
- `sbatch/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS.sbatch`
- `sbatch/M21-GENERANNO-1P2B-CRF-SCREEN-s1-shared.sbatch`
- `sbatch/M12C-GENERANNO-FAIR-CHALLENGER-SMOKE.sbatch`
- `sbatch/M14-ANIMAL-DISTANCE-NEGCTRL.sbatch`
- `sbatch/M12B-SAMEPANEL-BASELINES-EXTERNAL.sbatch`
- `sbatch/M13-DISTANCE-GENERALIZATION-SCAN.sbatch`
- `scripts/run_TA-COHERENCE-FIX-M5.sbatch`
- `scripts/run_BASE-TIBERIUS-PILOT-M1.sbatch`
- `scripts/run_screen_ref_M7.sbatch`
- `scripts/run_TA-FOUNDATION-DECODER-M4.sbatch`
- `scripts/run_FP-SEGMENTNT-PROBE-M1.sbatch`
- `scripts/run_BASE-HELIXER-SAC-DMEL-SMOKE-M1.sbatch`
- `scripts/run_BASE-HELIXER-SAC-SMOKE-M1.sbatch`
- `scripts/run_BASE-ANNEVO-SAC-DMEL-SMOKE-M1.sbatch`
- `scripts/run_screen_ref.sbatch`
- `ARCHITECTURE.md` (framework-only architecture guide)
