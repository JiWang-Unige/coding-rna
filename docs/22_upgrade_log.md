# Framework Upgrade Log / 框架升级记录

> 由 `/framework-upgrade` 维护。记录 auto-research 框架自身的版本迁移、兼容修复、新 skill/doc/hook 引入，以及哪些研究内容被保留。

## Upgrade Entries

### Upgrade 2026-06-14 — coding-rna v3.5+ → v4.1
- Reason: 用户请求将 `/home/users/j/jwang/auto-research-portable-v4.1.tar.gz` 的框架升级到 `/home/users/j/jwang/coding-rna`，以获得 v4.1 skills/hooks/docs、route reset、code-review gate、evaluator/baseline/code-review/upgrade 中央账本。
- Files changed: refreshed framework scripts/hooks, `.claude/`, `.codex/`, `.agents/`, `agents/`, `README.auto-research.md`, `ARCHITECTURE.md`, `PROJECT_STRUCTURE.md`, `cluster_config.yaml.example`, `.mcp.json.example`, `RUN_PROMPT.codex.md`; merged `CLAUDE.md` with v4.1 generic sections while preserving project-specific §0-2 and §12-14; regenerated `AGENTS.md`.
- Skills added/updated: v4.1 skill layer from the tarball, including `framework-upgrade`, `route-reset`, `code-review-gate`, `artifact-registry`, `publication-plan`, `pipeline-blueprint`, `sota-randomized`, `workspace-matrix`, and updated B-stage skills.
- Docs added/updated: seeded v4.1 docs `11-22`; populated `docs/11_master_plan.md`; migrated legacy `docs/11_evaluator_contract.md` into `docs/19_evaluator_contract.md`; migrated legacy `docs/12_baseline_reproduction.md` into `docs/20_baseline_reproduction.md`; appended this upgrade log and evidence register entry.
- Scripts/hooks changed: refreshed shipped v4.1 scripts without deleting project-specific scripts in `scripts/`; added/updated `pre_submit_gate.py`, `research_flow_guard.py`, `validate_codex_skills.py`, `validate_stage_c.py`, `artifact_registry.py`, `workspace_matrix.py`, `next_evidence_id.py`, `note_gate.py`, `sota_seed_matrix.py`, and hook scripts including `precompact_flush.sh` and `stage_flow_nudge.sh`.
- Research content preserved: `ACTIVE_GOAL.json`, `cluster_config.yaml`, `secrets.env`, `.mcp.json`, `docs/00-10`, old `docs/11_evaluator_contract.md`, old `docs/12_baseline_reproduction.md`, `docs/experiments/`, `goals/`, `refs/`, `wiki/`, `data/`, `outputs/`, `runs/`, `logs/`, `software_outputs/`, project-specific scripts and sbatch files were not deleted.
- Compatibility checks: installer completed `build_codex_skills.py`, `validate_codex_skills.py`, and `sync_agents_md.sh`; post-upgrade `python3 scripts/context_pack.py --purpose iterate` passed; `python3 scripts/validate_codex_skills.py .` passed; `python3 scripts/research_flow_guard.py . --format markdown` ran and returned expected research-state blockers rather than installation errors.
- Required follow-up: `research_flow_guard` recommends `/configure-project` because `ACTIVE_GOAL.json` remains `status=draft` and canonical training/eval templates still contain placeholders; it also reports evaluator/baseline readiness blockers because v4.1 structured checks require fully filled `docs/19/20` tables or reproduction ledger rows. Consider `$artifact-registry` audit because the project has many historical outputs and project-specific scripts.
- Rollback note: no git repository is present; installer backups use `.backup-20260614-214655`, and pre-merge driver backup is `CLAUDE.md.framework-upgrade-backup-20260614-214748`. Restore those files/directories if the v4.1 layer needs to be backed out.

### Upgrade 2026-06-13 — v4.0 → v4.1
- Reason: 补齐框架升级、同项目重开线、代码审前闸、evaluator/baseline 中央留档，以及人闸前通俗摘要规则。
- Files changed: `CLAUDE.md`, `README.auto-research.md`, `ARCHITECTURE.md`, `.claude/skills/*`, `scripts/*`, `docs/15-22`, install/sync scripts, hooks.
- Skills added/updated: added `/framework-upgrade`, `/route-reset`, `/code-review-gate`; updated `/implement`, `/pursue`, `/goal-prompt`, `/reproduce-baselines`, `/benchmark-roadmap`, `/artifact-registry`, `/master-plan`, `/council`, `/reframe`, `/workspace-matrix`, `/note-gate`.
- Docs added/updated: added `docs/19_evaluator_contract.md`, `docs/20_baseline_reproduction.md`, `docs/21_code_review_log.md`, `docs/22_upgrade_log.md`; updated `docs/15/16/17`.
- Scripts/hooks changed: `context_pack.py`, `research_flow_guard.py`, `iter_ledger.py`, `artifact_registry.py`, `note_gate.py`, `guard_paths.sh`, `submit_guard.sh`, `iter_record_nudge.sh`, `research_bootstrap.sh`, `install.sh`, `sync_agents_md.sh`.
- Research content preserved: no research result, refs, wiki, runs, outputs, data, secrets, ACTIVE_GOAL, or cluster_config was overwritten.
- Compatibility checks: `build_codex_skills.py`, `validate_codex_skills.py`, `sync_agents_md.sh`, `context_pack.py --purpose iterate`, `research_flow_guard.py`, `iter_ledger.py`, and new skill quick validators completed; false-positive stage guards were fixed.
- Required follow-up: review `git status`, then make one lightweight initial commit if the tracked file set looks right. No automatic commit was made.
- Rollback note: revert this upgrade by restoring the touched framework files from git/backup once git is initialized or from external backup; do not delete docs/19-22 if they already contain later project evidence.

### Upgrade <date> — <from> → <to>
- Reason:
- Files changed:
- Skills added/updated:
- Docs added/updated:
- Scripts/hooks changed:
- Research content preserved:
- Compatibility checks:
- Required follow-up:
- Rollback note:

### Upgrade 2026-07-01 — v4.1 → v4.2 lane-split
- Reason: User requested the same v4.2 lane-split architecture refresh as the TE project, using `/home/users/j/jwang/auto-research-portable-v4.2-lane-split-20260630.tar.gz`, to add the middle-lane workflow layer without overwriting coding-rna research state.
- Files changed: refreshed framework entry docs/templates `README.auto-research.md`, `README.md`, `ARCHITECTURE.md`, `PROJECT_STRUCTURE.md`, `RUN_PROMPT.codex.md`, `.gitignore`, `.mcp.json.example`, `cluster_config.yaml.example`, `secrets.env.example`, `install.sh`; refreshed `.claude/`, `.codex/`, `.agents/`, `agents/openai.yaml`, and framework scripts/hooks. `CLAUDE.md` was merged rather than overwritten: v4.2 framework sections were adopted while preserving coding-rna project-specific §0-2 and §12-14; `AGENTS.md` was regenerated from the merged `CLAUDE.md`.
- Skills added/updated: added `$evidence-sprint`, `$capability-pursue`, and `$review-board`; regenerated `.agents/skills` and `.codex/skills` from canonical `.claude/skills`.
- Docs added/updated: seed-if-absent created `docs/23_review_board.md` and `docs/24_sprint_pursue_ledger.md`; `docs/00-22*.md` research content was preserved except for this upgrade-log append and the evidence-register row documenting the upgrade.
- Scripts/hooks changed: refreshed v4.2 framework scripts from the release tarball, including `context_pack.py`, `research_flow_guard.py`, `validate_codex_skills.py`, `validate_stage_c.py`, and hook scripts; preserved project-specific experiment/training scripts already present in `scripts/`.
- Research content preserved: `ACTIVE_GOAL.json`, `cluster_config.yaml`, `secrets.env`, `.mcp.json`, `docs/00-22`, legacy `docs/11_evaluator_contract.md` and `docs/12_baseline_reproduction.md`, `docs/experiments/`, `refs/`, `wiki/`, `data/`, `configs/`, `outputs/`, `runs/`, `reports/`, `logs/`, `software_outputs/`, project-specific scripts, and sbatch files were not deleted or reset.
- Compatibility checks: `python3 scripts/build_codex_skills.py .` completed with total cross-agent description budget 6157 chars; `python3 scripts/validate_codex_skills.py .` passed for `.agents/skills` and `.codex/skills`; `bash scripts/sync_agents_md.sh` regenerated `AGENTS.md`; `python3 scripts/context_pack.py --purpose iterate` completed; final `validate_codex_skills.py` passed with 37 skills in each Codex-visible layer. `python3 scripts/research_flow_guard.py . --format markdown` returned the pre-existing research-state warning that the project is not fully configured/evaluator-ready and recommends `$configure-project`; this is not an upgrade validation failure.
- Required follow-up: no git repository is present in this workspace, so no commit/status review was performed. Consider `$configure-project` later if the project contract should be re-materialized for v4.2 guard semantics; use `$evidence-sprint` for 1-2 step evidence checks, `$capability-pursue` for bounded 2-5 round component work, and `$review-board` for non-result independent design review.
- Rollback note: restore framework files from existing timestamped backups or the previous v4.1 source if needed; do not delete `docs/23_review_board.md` or `docs/24_sprint_pursue_ledger.md` after they start accumulating project evidence.

## Compatibility Decisions
| Date | Decision | Reason | Affected files | Revisit condition |
|---|---|---|---|---|
| 2026-06-13 | Initialize git, but do not commit automatically | Framework/docs/scripts are now complex enough to need version history; heavy artifacts, secrets, cloned repos, PDFs, and runtime outputs remain ignored. | `.git/`, `.gitignore`, `.git/info/exclude`, `docs/17_parallel_workspace.md` | Revisit before first worktree cohort or if `.gitignore` would include sensitive/heavy files. |
| 2026-06-13 | Put per-run generated scripts under `scripts/experiments/<exp_id>/` | The existing contract separated reusable scripts and sbatch/results, but did not name a home for one-off wrappers that affect a run. | `docs/16_artifact_registry.md`, `PROJECT_STRUCTURE.md`, `scripts/artifact_registry.py`, `/artifact-registry`, `/implement`, `CLAUDE.md`, `AGENTS.md`, `README.auto-research.md`, `ARCHITECTURE.md` | Revisit if pipeline stages need a different namespace; pipeline DAG code still belongs under `pipelines/<pipeline_id>/`. |
| 2026-06-14 | Accept metric summaries from `reports/<exp_id>.json` or `outputs/<exp_id>/metrics/metrics.json` | The v4.1 template prefers `reports/`, but `coding-rna` already has many valid metrics under `outputs/<exp_id>/metrics/` and `CLAUDE.md §12` names that as the canonical project path. | `scripts/artifact_registry.py`, `docs/16_artifact_registry.md`, `PROJECT_STRUCTURE.md`, `docs/05_todo.md` | Revisit only if future training code standardizes on `reports/` and a migration is explicitly desired. |
| 2026-06-14 | Current submission mode is `on_cluster` on Baobab login node | Live probe shows `sinfo`/`sbatch`/`sacct`/`conda` in PATH on `login1.baobab`; old `remote_ssh`/`ssh baobab` contract was for running the driver off-cluster. | `CLAUDE.md`, `AGENTS.md`, `cluster_config.yaml`, `ACTIVE_GOAL.json`, `docs/11_master_plan.md` | Revisit if the driver runs from a non-cluster workstation again. |
