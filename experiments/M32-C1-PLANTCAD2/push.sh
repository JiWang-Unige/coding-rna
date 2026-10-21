#!/bin/bash
# push local M32-C1-PLANTCAD2 (scripts/tests; not results) to ont ~/coding-rna-reorg/experiments/M32-C1-PLANTCAD2
cd /Users/jiwang/Desktop/coding-rna-reorg/experiments || exit 1
COPYFILE_DISABLE=1 tar c --exclude results --exclude __pycache__ M32-C1-PLANTCAD2 | ssh baobab 'ssh lab "mkdir -p ~/coding-rna-reorg/experiments && tar x -C ~/coding-rna-reorg/experiments"'
