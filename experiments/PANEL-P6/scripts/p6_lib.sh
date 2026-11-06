# sourced by the PANEL-P6 GPU drivers (p6_native.sh, p6_carbona.sh, p6_vipsania.sh, p6_queue.sh)
# Paths, the event log, and the GPU gate.  Gate (checked before every GPU step, i.e. before every genome of every tool):
#   1. $P6/GO exists (P6_GO overrides the file, e.g. GO_SMOKE for the smoke tests; written by hand once the ont slot is agreed with the maintenance session; $P6/STOP stops the queue between steps);
#   2. no TE GPU claim file ($TE_CLAIM, written by the TE session before its GPU work and removed after it) and
#      no GPU process of the TE project (user rule 2026-11-07: P6 runs only when no TE job runs; never pre-empt): a GPU process is TE when its cwd or args contain
#      /work/te, ab-initio-TE, TE-GLM or te_benchmark, or when it runs inside an Apptainer/Singularity container that is not ours (our containers carry PANEL-P6 in their args);
#      only nvidia-smi compute processes are examined (TE's CPU-only Apptainer jobs such as BRAKER never block);
#   3. no other coding-rna GPU job and no live M34 / M39 queue driver (one coding-rna GPU job at a time; P6 does not slip into the gap between two queued arms)
#      unless P6_SHARE_WITH_CODING_RNA=1;
#   4. at least NEED MiB of GPU memory free (2 consecutive checks, 30 s apart).
# While waiting the gate logs one gate_wait event per hour.
P6=${P6_OUT:-$HOME/coding-rna/outputs/PANEL-P6}
P6S=$HOME/coding-rna-reorg/experiments/PANEL-P6/scripts
A=$P6/assets
G=$P6/genomes
LOGS=${P6_LOGS:-$HOME/coding-rna/logs}
EV=$P6/events.jsonl
APPT=$HOME/envs/apptainer/bin/apptainer
PYG=$HOME/envs/generanno/bin/python
CPUS=${P6_CPUS:-2,4,6,8,10,12,14,16-24}               # the cores of the coding-rna GPU steps (M34 / M38)
PCPUS=${P6_POST_CPUS:-25-31}                          # CPU-only side work (Carbon-A post-processing, conversion)
GUARD="python3.12 -I $HOME/coding-rna-reorg/experiments/M38-VIPSANIA/scripts/vip_guard.py"
EVT="python3.12 -I $HOME/coding-rna-reorg/experiments/M38-VIPSANIA/scripts/vip_event.py --file $EV"
KEYS_BY_SIZE="dryas_octopetala juncus_effusus acorus_gramineus eschscholzia_californica rhynchospora_breviuscula zizania_latifolia akebia_trifoliata amborella_trichopoda hakea_prostrata themeda_triandra"
mkdir -p $P6 $LOGS
ev() { $EVT "$@"; }
gpu_free() { nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -dc 0-9; }
fsize() { stat -c %s "$1" 2>/dev/null || echo 0; }
disk_free_gb() { df -BG --output=avail $HOME | tail -1 | tr -dc 0-9; }

gpu_blockers() {   # prints one line per GPU process that blocks us (TE, other coding-rna jobs); empty = none
  local pid cwd args par pargs
  for pid in $(nvidia-smi --query-compute-apps=pid --format=csv,noheader,nounits 2>/dev/null | tr -dc '0-9\n'); do
    cwd=$(readlink /proc/$pid/cwd 2>/dev/null); args=$(tr '\0' ' ' < /proc/$pid/cmdline 2>/dev/null)
    case "$cwd $args" in *PANEL-P6*) continue;; esac                      # ours
    case "$cwd $args" in */work/te*|*ab-initio-TE*|*TE-GLM*|*te_benchmark*) echo "TE pid=$pid $args" | cut -c1-200; continue;; esac
    if [ -z "${P6_SHARE_WITH_CODING_RNA:-}" ]; then                         # set only if the maintenance session agrees to share the card with M34 / M39
      case "$cwd $args" in *M34-NTV2-SCALEUP*|*M39-INSECT-C2*|*m34_*|*m39_*) echo "coding-rna pid=$pid $args" | cut -c1-200; continue;; esac
    fi
    par=$pid                                                                # a container that is not ours counts as TE (TE runs Tiberius / BRAKER in Apptainer)
    for _ in 1 2 3 4 5 6; do
      par=$(awk '/^PPid:/{print $2}' /proc/$par/status 2>/dev/null); [ -z "$par" ] || [ "$par" -le 1 ] && break
      pargs=$(tr '\0' ' ' < /proc/$par/cmdline 2>/dev/null)
      case "$pargs" in *PANEL-P6*) break;; *apptainer*|*singularity*|*starter*) echo "container pid=$pid parent=$par $pargs" | cut -c1-200; break;; esac
    done
  done
  if [ -z "${P6_SHARE_WITH_CODING_RNA:-}" ]; then                           # the coding-rna queue drivers themselves (no GPU process between two arms, but the queue goes on)
    pgrep -af 'm34_[a]fter_|m34_[s]tage.sh|m34_[s]upervise.sh|m39_[a]fter_|m39_[r]un_all' | cut -c1-160 | sed 's/^/coding-rna queue /'
  fi
}

TE_CLAIM=${P6_TE_CLAIM:-$HOME/ab-initio-TE/outputs/TE_GPU_CLAIM}   # TE's GPU claim (TE session, 2026-11-07): present while TE owns the ont GPU

gate() {   # NEED_MiB STEP : blocks until the gate is open; returns 1 if $P6/STOP appears
  local need=$1 step=$2 ok=0 last=0 f b
  while :; do
    [ -e $P6/STOP ] && { ev event=stopped step=$step; return 1; }
    b=$(gpu_blockers); f=$(gpu_free); f=${f:-0}
    [ -e $TE_CLAIM ] && b="${b}TE_GPU_CLAIM: $(head -c 300 $TE_CLAIM | tr '\n' ' ')"$'\n'
    if [ -e ${P6_GO:-$P6/GO} ] && [ -z "$b" ] && [ "$f" -ge "$need" ] && [ "$(disk_free_gb)" -ge ${P6_MIN_DISK_GB:-120} ]; then
      ok=$((ok+1)); [ $ok -ge 2 ] && return 0; sleep 30; continue
    fi
    ok=0
    if [ $(( $(date +%s) - last )) -ge 3600 ]; then
      last=$(date +%s)
      ev event=gate_wait step=$step go=$([ -e ${P6_GO:-$P6/GO} ] && echo 1 || echo 0) gpu_free_mb=$f need_mb=$need disk_free_gb=$(disk_free_gb) "blockers=$(echo "$b" | tr '\n' '|' | cut -c1-400)"
    fi
    sleep 120
  done
}

run_step() {   # NAME TIMEOUT_S MAX_RSS_GB CMD... : one command under the RAM / timeout / GPU-memory guard; logs ~/coding-rna/logs/p6-<name>.{out,err}; resources $P6/steps/<name>.res.json
  local NAME=$1 TO=$2 RSS=$3; shift 3
  mkdir -p $P6/steps
  local OUT=$LOGS/p6-$NAME.out ERR=$LOGS/p6-$NAME.err RES=$P6/steps/$NAME.res.json
  rm -f $RES
  ev event=step_start step=$NAME "gpu_free_mb=$(gpu_free)" "cmd=$*"
  $GUARD --gpu --max-rss-gb $RSS --timeout-s $TO --interval 15 --out $RES -- nice -n 10 taskset -c $CPUS "$@" >> $OUT 2>> $ERR < /dev/null
  local rc=$?
  $EVT --merge $RES event=step_end step=$NAME
  return $rc
}
