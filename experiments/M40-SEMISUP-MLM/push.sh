#!/bin/bash
# push local M40-SEMISUP-MLM (scripts/tests; not results) to ont ~/coding-rna-reorg/experiments/M40-SEMISUP-MLM
cd /Users/jiwang/Desktop/coding-rna-reorg/experiments || exit 1
COPYFILE_DISABLE=1 tar c --exclude results --exclude __pycache__ M40-SEMISUP-MLM | ssh baobab 'ssh lab "mkdir -p ~/coding-rna-reorg/experiments && tar x -C ~/coding-rna-reorg/experiments"'
