#!/bin/bash
# run on ont: bash tests/run_tests.sh   (generanno python; CPU only, cores 25-31, nice)
cd "$(dirname "$0")"
PY=$HOME/envs/generanno/bin/python
for t in test_map_roundtrip.py test_projection_flip.py test_lambda0_identity.py; do
  echo "== $t"; nice -n 10 taskset -c 25-31 $PY $t || exit 1
done
echo ALL_TESTS_PASSED
