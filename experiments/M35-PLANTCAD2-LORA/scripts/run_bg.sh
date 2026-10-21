#!/bin/bash
# usage: run_bg.sh <name> <command...>   detached (setsid nohup) with all fds redirected; log ~/coding-rna/logs/<name>-<ts>.{out,err}; prints "PID LOG"
NAME=$1; shift
LOG=$HOME/coding-rna/logs/$NAME-$(date +%Y%m%d-%H%M%S)
setsid nohup "$@" > $LOG.out 2> $LOG.err < /dev/null &
echo "$! $LOG"
