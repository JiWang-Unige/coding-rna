#!/bin/bash
# Control-flow test of vip_devb.sh with a stub `vipsania` (no GPU, no TensorFlow, no model): events, sentinels, resumability, refusal without READY.
# Runs under a throw-away root (M38_ROOT), never touches the real $HOME/coding-rna/outputs/M38-VIPSANIA/READY.
#   bash test_devb_stub.sh            (on ont; python3.12 + nvidia-smi query only)
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SC=$HERE/../scripts
T=${1:-$HOME/coding-rna/outputs/M38-VIPSANIA/_stub_test}
rm -rf $T; mkdir -p $T/genomes $T/vipsania_cache/models/j9m0cdmk $T/logs $T/bin
python3.12 -I - <<PY
import random
random.seed(1)
def fa(path, contigs):
    with open(path, "w") as h:
        for name, n in contigs:
            h.write(f">{name} test\n"); s = "".join(random.choice("ACGT") for _ in range(n))
            for i in range(0, n, 80): h.write(s[i:i+80] + "\n")
import json, os
man = {}
for g, c in (("arabidopsis_thaliana", [("NC_003070.9", 7100000)]), ("dryas_octopetala", [("NC_141367.1", 5000)]), ("oryza_sativa", [("NC_089041.1", 6000)]), ("eschscholzia_californica", [("NW_028238309.1", 7000)])):
    p = f"$T/genomes/{g}.nuclear.fa"; fa(p, c); man[g] = dict(bytes=os.path.getsize(p))
json.dump(man, open("$T/genomes/genomes_manifest.json", "w"))
open("$T/vipsania_cache/models/j9m0cdmk/config.json", "w").write("{}"); open("$T/vipsania_cache/models/j9m0cdmk/latest_checkpoint.weights.h5", "w").write("w")
PY
(cd $T/vipsania_cache/models/j9m0cdmk && sha256sum config.json latest_checkpoint.weights.h5 > $T/model_j9m0cdmk.sha256)
cat > $T/bin/vipsania <<'STUB'
#!/bin/bash
# stub: vipsania annotate MODEL FASTA -o OUT [--finetune ...]
out=""; ft=0; while [ $# -gt 0 ]; do case "$1" in -o) out=$2; shift;; --finetune) ft=1;; esac; shift; done
echo "Using 1 GPU: stub"
mkdir -p $(dirname $out)
for i in 1 2 3; do printf 'ctg\tVipsania\tgene\t%d\t%d\t.\t+\t.\tID=g%d\n' $i $((i+300)) $i; done > $out
echo "> time: 0.01 minutes" > ${out%.*}.log
[ "${STUB_FAIL_FT:-0}" = 1 ] && [ $ft = 1 ] && { echo "ResourceExhaustedError: OOM when allocating tensor" >&2; exit 1; }
exit 0
STUB
chmod +x $T/bin/vipsania
export VIP_GPU_FREE_MB=24000 M38_ROOT=$T VIP_BIN=$T/bin/vipsania VIP_LOGS=$T/logs VIP_SKIP_SCORE=1 VIP_MIN_GENES=1 VIP_SMOKE_MIN_GENES=1 VIP_CPUS=0-3 VIP_FT_TIMEOUT_H=1 VIP_ANN_TIMEOUT_H=1
ok=1
chk() { if eval "$2"; then echo "[OK ] $1"; else echo "[FAIL] $1"; ok=0; fi; }
bash $SC/vip_devb.sh; rc=$?
chk "refused without READY (rc 11)" "[ $rc = 11 ] && grep -q no_READY $T/events.jsonl"
touch $T/READY
bash $SC/vip_devb.sh; rc=$?
chk "run completes (rc 0)" "[ $rc = 0 ]"
chk "DEVB_DONE written" "[ -s $T/DEVB_DONE ]"
chk "per-genome sentinels" "[ -e $T/runs/arabidopsis_thaliana/DONE ] && [ -e $T/runs/eschscholzia_californica/DONE ]"
chk "step sentinels (pre + ft per genome)" "[ \$(ls $T/runs/*/pre/pre_DONE $T/runs/*/ft/ft_DONE | wc -l) = 8 ]"
chk "genome order smallest first" "[ \"\$(grep -o 'event\": \"genome_done\", \"genome\": \"[a-z_]*' $T/events.jsonl | sed 's/.*: \"//' | tr '\\n' ' ')\" = 'arabidopsis_thaliana dryas_octopetala oryza_sativa eschscholzia_californica ' ]"
chk "GPU smoke recorded finetune_B and projection" "grep -q gpu_smoke_done $T/events.jsonl && [ -s $T/ft_B.txt ]"
chk "run.json has seconds / model / genes" "python3.12 -I -c \"import json; r=json.load(open('$T/runs/oryza_sativa/ft/run.json')); assert r['model']=='j9m0cdmk' and r['genes']==3 and r['seconds']>=0\""
chk "logs named m38-vipsania-*.out/err" "ls $T/logs/m38-vipsania-arabidopsis_thaliana_ft.out $T/logs/m38-vipsania-gpusmoke_pre.err >/dev/null"
n=$(wc -l < $T/events.jsonl)
bash $SC/vip_devb.sh; rc=$?
chk "second invocation exits at once (DEVB_DONE)" "[ $rc = 0 ] && grep -q already_done $T/events.jsonl"
# resumability: drop DEVB_DONE and the last genome's ft sentinels -> only that part is redone
rm -f $T/DEVB_DONE $T/runs/eschscholzia_californica/DONE $T/runs/eschscholzia_californica/ft/ft_DONE
bash $SC/vip_devb.sh; rc=$?
chk "resume: finished genomes skipped, last ft redone" "[ $rc = 0 ] && grep -q '\"event\": \"genome_skipped\", \"genome\": \"oryza_sativa\"' $T/events.jsonl && [ -e $T/runs/eschscholzia_californica/ft/ft_DONE ] && [ -s $T/DEVB_DONE ]"
# projection guard: an (artificially tiny) limit stops the job after the smoke, before any genome
rm -rf $T/runs $T/DEVB_DONE $T/gpu_smoke
VIP_MAX_PROJ_H=0.001 bash $SC/vip_devb.sh; rc=$?
chk "projection above the limit aborts with rc 25 before any genome" "[ $rc = 25 ] && [ ! -e $T/runs/arabidopsis_thaliana/DONE ] && grep -q 'stage\": \"projection' $T/events.jsonl"
# OOM ladder: finetune OOMs for B=2 -> B=1 recorded (stub fails every ft run, so the script must give up cleanly with a nonzero rc)
rm -rf $T/runs/dryas_octopetala $T/DEVB_DONE $T/gpu_smoke; echo 2 > $T/ft_B.txt
STUB_FAIL_FT=1 bash $SC/vip_devb.sh; rc=$?
chk "persistent finetune failure aborts nonzero (no DEVB_DONE)" "[ $rc != 0 ] && [ ! -e $T/DEVB_DONE ] && grep -q devb_aborted $T/events.jsonl"
[ $ok = 1 ] && echo "STUB TEST PASSED" || { echo "STUB TEST FAILED"; exit 1; }
