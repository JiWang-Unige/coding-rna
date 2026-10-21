#!/bin/bash
# push local M35-PLANTCAD2-LORA (scripts/tests; not results) to ont ~/coding-rna-reorg/experiments/M35-PLANTCAD2-LORA
cd /Users/jiwang/Desktop/coding-rna-reorg/experiments || exit 1
COPYFILE_DISABLE=1 tar c --exclude results --exclude __pycache__ M35-PLANTCAD2-LORA | ssh baobab 'ssh lab "mkdir -p ~/coding-rna-reorg/experiments && tar x -C ~/coding-rna-reorg/experiments"'
