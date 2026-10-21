#!/bin/bash
# push local M34-NTV2-SCALEUP (scripts/tests; not results) to ont ~/coding-rna-reorg/experiments/M34-NTV2-SCALEUP
cd /Users/jiwang/Desktop/coding-rna-reorg/experiments || exit 1
COPYFILE_DISABLE=1 tar c --exclude results --exclude __pycache__ M34-NTV2-SCALEUP | ssh baobab 'ssh lab "mkdir -p ~/coding-rna-reorg/experiments && tar x -C ~/coding-rna-reorg/experiments"'
