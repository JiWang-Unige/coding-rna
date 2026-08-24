# Pipeline Blueprint / 流程化研究推进图

> 由 `/pipeline-blueprint` 维护。用于“已有既定思路、raw data 或分析流程，需要把 pipeline 建稳并产出可投稿证据”的阶段。它不是 blind iteration，而是 DAG 化推进。

## 0. Pipeline identity
- Pipeline name:
- Purpose:
- Input raw data:
- Final outputs:
- Claim supported:

## 1. DAG overview

```text
raw_data → QC → preprocessing → feature/model/statistical analysis → validation → figures/tables
```

| Stage | Purpose | Input | Output | Script/software | Parameters/config | QC gate | Status |
|---|---|---|---|---|---|---|---|
| S1 |  |  |  |  |  |  | TODO |

## 2. IO contracts

| Artifact | Path pattern | Producer | Consumer | Required metadata/hash | Retention policy |
|---|---|---|---|---|---|
| raw data | `data/raw/...` |  |  |  | keep |
| processed data | `data/processed/...` |  |  |  | keep |
| external software output | `software_outputs/<tool>/<run_id>/...` |  |  |  | keep summary + logs |

## 3. External software calls

| Tool | Version/container/env | Command template | Output dir | Log file | Failure handling |
|---|---|---|---|---|---|
|  |  |  |  |  |  |

## 4. Validation and sensitivity plan

| Check | Why | Method | Pass criterion | Status |
|---|---|---|---|---|
|  |  |  |  | TODO |

## 5. Pipeline execution ledger

| Run ID | Date | Stage(s) | Input hash | Config | Output | Status | Notes |
|---|---|---|---|---|---|---|---|

## 6. Handoff
- Next stage to execute:
- Blocking input needed:
- User decision needed:
