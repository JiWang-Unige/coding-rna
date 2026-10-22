#!/bin/bash
# push local M37-JOINT-DECODE-SMOKE (scripts/tests; not results) to ont ~/coding-rna-reorg/experiments/M37-JOINT-DECODE-SMOKE
cd /Users/jiwang/Desktop/coding-rna-reorg/experiments || exit 1
COPYFILE_DISABLE=1 tar c --exclude results --exclude __pycache__ M37-JOINT-DECODE-SMOKE | ssh baobab 'ssh lab "mkdir -p ~/coding-rna-reorg/experiments && tar x -C ~/coding-rna-reorg/experiments"'
