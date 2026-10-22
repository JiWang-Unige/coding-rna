#!/bin/bash
# push local M36-BOOTSTRAP-ADAPT (scripts/tests; not results) to ont ~/coding-rna-reorg/experiments/M36-BOOTSTRAP-ADAPT
cd /Users/jiwang/Desktop/coding-rna-reorg/experiments || exit 1
COPYFILE_DISABLE=1 tar c --exclude results --exclude __pycache__ M36-BOOTSTRAP-ADAPT | ssh baobab 'ssh lab "mkdir -p ~/coding-rna-reorg/experiments && tar x -C ~/coding-rna-reorg/experiments"'
