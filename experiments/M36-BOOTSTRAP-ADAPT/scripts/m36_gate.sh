#!/bin/bash
# Gate functions of the M36 queue (sourced; also tested).  Same logic as the P2b driver (m34_p2b.sh): the ab-initio-TE project shares this GPU under the same account, so only
# GPU compute processes that belong to coding-rna (command line or working directory contains 'coding-rna') block.  The pgrep patterns contain a bracket so that pgrep never
# matches its own command line.   Variables: M36_ME (user), P2B_DONE (sentinel path), M36_GATE_SKIP_PROCS=1 (dry runs: skip the process / GPU scan)
M36_ME=${M36_ME:-$(id -un)}

m36_gpu_blockers() {     # prints the pids of coding-rna GPU compute processes of this user (empty = none)
  local g="" p
  for p in $(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null); do
    [ "$(ps -o user= -p $p 2>/dev/null | tr -d ' ')" = "$M36_ME" ] || continue
    if tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | grep -q 'coding-rna' || readlink /proc/$p/cwd 2>/dev/null | grep -q 'coding-rna'; then g="$g $p"; fi
  done
  echo $g
}

m36_blockers() {         # reasons that block the START of the queue (empty = free): no P2B_DONE, other coding-rna jobs, coding-rna GPU processes
  local r=""
  [ -f "$P2B_DONE" ] || r="no P2B_DONE;"
  if [ "${M36_GATE_SKIP_PROCS:-0}" != 1 ]; then
    local pp; pp=$(pgrep -u $M36_ME -f 'm34_trai[n]|m34_infe[r]|m34_p2[b]|m34_stag[e]|m34_queu[e]|m34_supervis[e]|m34_smok[e]' | tr '\n' ' ')
    [ -n "$pp" ] && r="$r processes: $pp;"
    local g; g=$(m36_gpu_blockers)
    [ -n "$g" ] && r="$r gpu pids: $g;"
  fi
  echo "$r" | tr -d '"\\'
}
