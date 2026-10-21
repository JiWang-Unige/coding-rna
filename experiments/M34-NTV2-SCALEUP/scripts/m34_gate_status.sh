#!/bin/bash
# one-line status for the waiter: gate state, M32 progress
ls ~/coding-rna/outputs/M32-C1-PLANTCAD2/M32_DONE ~/coding-rna/outputs/M33-NTV2-FOLLOWUP/M33_DONE 2>&1 | sed 's/^.*No such file.*$/missing/' | tr '\n' ' '; echo
tail -n 1 ~/coding-rna/outputs/M32-C1-PLANTCAD2/runA/training.jsonl 2>/dev/null | cut -c1-80
