#!/bin/bash
# push local PANEL-PREP (scripts/tests; not results) to ont ~/coding-rna-reorg/experiments/PANEL-PREP
cd /Users/jiwang/Desktop/coding-rna-reorg/experiments || exit 1
COPYFILE_DISABLE=1 tar c --exclude results --exclude __pycache__ --exclude .pytest_cache PANEL-PREP | ssh baobab 'ssh lab "mkdir -p ~/coding-rna-reorg/experiments && tar x -C ~/coding-rna-reorg/experiments"'
