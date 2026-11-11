#!/bin/bash
# push local M41-INTRON-LENGTH (scripts/tests/README; not results) to ont ~/coding-rna-reorg/experiments/M41-INTRON-LENGTH
# (route: ont is reachable only via the baobab jump; do NOT add </dev/null to either ssh, the tar stream would be lost)
cd /Users/jiwang/Desktop/coding-rna-reorg/experiments || exit 1
COPYFILE_DISABLE=1 tar c --exclude results --exclude __pycache__ M41-INTRON-LENGTH | ssh baobab 'ssh lab "mkdir -p ~/coding-rna-reorg/experiments && tar x -C ~/coding-rna-reorg/experiments"'
