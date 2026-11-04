#!/bin/bash
# push local M38-VIPSANIA (scripts/README; not results) to ont ~/coding-rna-reorg/experiments/M38-VIPSANIA
# (route: ont is reachable only via the baobab jump; do NOT add </dev/null to either ssh, the tar stream would be lost)
cd /Users/jiwang/Desktop/coding-rna-reorg/experiments || exit 1
COPYFILE_DISABLE=1 tar c --exclude results --exclude __pycache__ M38-VIPSANIA | ssh baobab 'ssh lab "mkdir -p ~/coding-rna-reorg/experiments && tar x -C ~/coding-rna-reorg/experiments"'
