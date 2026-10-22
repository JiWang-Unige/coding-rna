#!/bin/bash
# Independent-panel download (public NCBI FTP; CPU/network only; same logic as M34 m34_download.sh).  Per assembly: genomic.fna.gz, genomic.gff.gz,
# protein.faa.gz, assembly_report.txt, assembly_stats.txt, md5checksums.txt, plus the NCBI datasets v2 dataset_report (annotation name/date/method, BUSCO).
# Verifies md5 against NCBI's md5checksums.txt; records local md5 + download date.   usage: panel_download.sh [species_key ...]
HERE=$(cd "$(dirname "$0")" && pwd)
RAW=$HOME/coding-rna/data/panel/raw
mkdir -p "$RAW"
KEYS=${@:-$(python3 -c "import json;print(' '.join(s['key'] for s in sorted(json.load(open('$HERE/panel_species.json')),key=lambda r:r['order'])))")}
for K in $KEYS; do
  read BASE ACC NAME < <(python3 - <<PY
import json
s=[x for x in json.load(open('$HERE/panel_species.json')) if x['key']=='$K'][0]
print(s['ftp_base'],s['accession'],s['assembly_name'])
PY
)
  D=$RAW/${ACC}_${NAME}; mkdir -p "$D"; DIR=$(dirname "$BASE"); STEM=$(basename "$BASE")
  echo "[$(date -Is)] $K -> $D"
  for SUF in _genomic.fna.gz _genomic.gff.gz _protein.faa.gz _assembly_report.txt _assembly_stats.txt; do
    [ -s "$D/$STEM$SUF" ] || nice -n 10 curl -fsSL --retry 5 --retry-delay 10 --max-time 3600 -o "$D/$STEM$SUF.part" "$BASE$SUF" && mv -f "$D/$STEM$SUF.part" "$D/$STEM$SUF" 2>/dev/null
  done
  [ -s "$D/md5checksums.txt" ] || curl -fsSL --retry 5 -o "$D/md5checksums.txt" "$DIR/md5checksums.txt"
  [ -s "$D/ncbi_dataset_report.json" ] || curl -fsSL --retry 5 -o "$D/ncbi_dataset_report.json" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$ACC/dataset_report"
  ( cd "$D"; : > local_md5.txt; ok=1
    for SUF in _genomic.fna.gz _genomic.gff.gz _protein.faa.gz _assembly_report.txt _assembly_stats.txt; do
      F=$STEM$SUF; L=$(md5sum "$F" | cut -d' ' -f1); R=$(grep " \./$F\$" md5checksums.txt | cut -d' ' -f1)
      echo "$L  $F  ncbi=$R  $([ "$L" = "$R" ] && echo OK || echo MISMATCH)" >> local_md5.txt; [ "$L" = "$R" ] || ok=0
    done
    date -Is > download_date.txt; echo "ftp_base=$BASE" >> download_date.txt
    [ $ok = 1 ] && echo "[$(date -Is)] $K md5 OK" || echo "[$(date -Is)] $K md5 MISMATCH"
  )
done
echo "[$(date -Is)] ALL_DONE"
