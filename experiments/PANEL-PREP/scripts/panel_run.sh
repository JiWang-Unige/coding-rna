#!/bin/bash
# ont driver (CPU, cores 25-31, nice): prep all 12 -> reference stats -> full-proteome homology audit -> result tables.  Logs to ~/coding-rna/outputs/PANEL-PREP/run.log
set -o pipefail
HERE=$(cd "$(dirname "$0")" && pwd); PY=$HOME/envs/generanno/bin/python; OUT=$HOME/coding-rna/outputs/PANEL-PREP; mkdir -p $OUT/results
cd $HERE
echo "[$(date -Is)] START prep"
$PY panel_prep.py prep || { echo "[$(date -Is)] PREP_FAILED"; exit 1; }
echo "[$(date -Is)] START stats"
$PY panel_prep.py stats || { echo "[$(date -Is)] STATS_FAILED"; exit 1; }
echo "[$(date -Is)] START homology"
$PY panel_homology.py 6 || { echo "[$(date -Is)] HOMOLOGY_FAILED"; exit 1; }
echo "[$(date -Is)] START tables"
$PY panel_tables.py $OUT/results && $PY panel_exposure.py $OUT/results/exposure_audit.json > $OUT/results/exposure_audit.txt
cp ~/coding-rna/data/panel/homology/homology_audit_full.json ~/coding-rna/data/panel/homology/homology_table.tsv $OUT/results/
echo "[$(date -Is)] PANEL_RUN_DONE"
