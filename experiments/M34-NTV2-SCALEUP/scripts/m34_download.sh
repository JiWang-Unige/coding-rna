#!/bin/bash
# M34 data download (public NCBI RefSeq FTP; CPU/network only).  Files per assembly: genomic.fna.gz, genomic.gff.gz, protein.faa.gz,
# assembly_report.txt, md5checksums.txt.  Verifies md5 against NCBI's md5checksums.txt and records local md5 + download date.
# usage: m34_download.sh [species_key ...]   (default: all in species.json)
HERE=$(cd "$(dirname "$0")" && pwd)
RAW=$HOME/coding-rna/data/m34_train/raw
mkdir -p "$RAW"
KEYS=${@:-$(python3 -c "import json;print(' '.join(s['key'] for s in json.load(open('$HERE/species.json'))))")}
for K in $KEYS; do
  read BASE ACC NAME < <(python3 - <<PY
import json
s=[x for x in json.load(open('$HERE/species.json')) if x['key']=='$K'][0]
print(s['ftp_base'],s['accession'],s['assembly_name'])
PY
)
  D=$RAW/${ACC}_${NAME}; mkdir -p "$D"; DIR=$(dirname "$BASE"); STEM=$(basename "$BASE")
  echo "[$(date -Is)] $K -> $D"
  for SUF in _genomic.fna.gz _genomic.gff.gz _protein.faa.gz _assembly_report.txt; do
    [ -s "$D/$STEM$SUF" ] || nice -n 10 curl -fsSL --retry 5 --retry-delay 10 --max-time 3600 -o "$D/$STEM$SUF.part" "$BASE$SUF" && mv "$D/$STEM$SUF.part" "$D/$STEM$SUF"
  done
  [ -s "$D/md5checksums.txt" ] || curl -fsSL --retry 5 -o "$D/md5checksums.txt" "$DIR/md5checksums.txt"
  ( cd "$D"; : > local_md5.txt; ok=1
    for SUF in _genomic.fna.gz _genomic.gff.gz _protein.faa.gz; do
      F=$STEM$SUF; L=$(md5sum "$F" | cut -d' ' -f1); R=$(grep " \./$F\$" md5checksums.txt | cut -d' ' -f1)
      echo "$L  $F  ncbi=$R  $([ "$L" = "$R" ] && echo OK || echo MISMATCH)" >> local_md5.txt; [ "$L" = "$R" ] || ok=0
    done
    date -Is > download_date.txt; echo "ftp_base=$BASE" >> download_date.txt
    [ $ok = 1 ] && echo "[$(date -Is)] $K md5 OK" || echo "[$(date -Is)] $K md5 MISMATCH"
  )
done
echo "[$(date -Is)] ALL_DONE"
