#!/bin/bash
# push local M33-NTV2-FOLLOWUP (scripts/tests/protocol; not results) to ont ~/coding-rna-reorg/experiments/M33-NTV2-FOLLOWUP
cd /Users/jiwang/Desktop/coding-rna-reorg/experiments || exit 1
COPYFILE_DISABLE=1 tar c --exclude results --exclude __pycache__ M33-NTV2-FOLLOWUP | ssh baobab 'ssh lab "mkdir -p ~/coding-rna-reorg/experiments && tar x -C ~/coding-rna-reorg/experiments"'
