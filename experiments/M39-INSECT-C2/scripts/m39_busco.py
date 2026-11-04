"""Protocol 29 §5.3 V7: uniform BUSCO (insecta_odb10, protein mode) of the annotated proteome of every M39 species, one longest protein per protein-coding gene
(the NCBI convention).  CPU only.  Step 1 writes data/m39_insect/busco/<key>.faa from the raw NCBI GFF3 (CDS protein_id -> mRNA -> gene) and protein.faa;
step 2 runs BUSCO 5.8.3 (pixi env ~/envs/busco) offline against the lineage downloaded with `busco --download insecta_odb10` (2024-01-08 release); results go to <OUT>/busco/<key>/ and a summary to <OUT>/busco/busco_summary.json.
usage: m39_busco.py faa [key ...] | run [cpus_per_job=4] [jobs=2] | summary
"""
import gzip
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import m39_common as C

PIXI = [str(Path.home() / '.pixi/bin/pixi'), 'run', '--manifest-path', str(Path.home() / 'envs/busco/pixi.toml')]   # BUSCO 5.8.3 (the env's script has a /usr/bin/env python3 shebang)
LINEAGE = 'insecta_odb10'
FAA = C.DATA / 'busco'
RES = C.OUT / 'busco'


def attrs(s):
    return dict(x.split('=', 1) for x in s.split(';') if '=' in x)


def longest_per_gene(row):
    raw = C.raw_dir(row)
    tx_gene, prot_tx, coding = {}, {}, set()
    with gzip.open(raw / f"{row['ftp_stem']}_genomic.gff.gz", 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            r = line.rstrip('\n').split('\t')
            if len(r) < 9:
                continue
            a = attrs(r[8])
            if r[2] == 'gene' and a.get('gene_biotype') == 'protein_coding':
                coding.add(a['ID'])
            elif r[2] == 'CDS' and 'protein_id' in a:
                prot_tx[a['protein_id']] = a['Parent'].split(',')[0]
            elif 'Parent' in a and 'ID' in a and r[2] not in ('exon', 'CDS'):
                tx_gene[a['ID']] = a['Parent'].split(',')[0]
    seqs, name, buf = {}, None, []
    with gzip.open(raw / f"{row['ftp_stem']}_protein.faa.gz", 'rt') as f:
        for line in f:
            if line[0] == '>':
                if name:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = ''.join(buf)
    best = {}
    for pid, s in seqs.items():
        g = tx_gene.get(prot_tx.get(pid))
        if g in coding and (g not in best or len(s) > len(seqs[best[g]])):
            best[g] = pid
    return {pid: seqs[pid] for pid in sorted(best.values())}, dict(proteins=len(seqs), coding_genes=len(coding), genes_with_protein=len(best))


def main():
    cmd = sys.argv[1]
    rows = C.roster()
    if cmd == 'faa':
        FAA.mkdir(parents=True, exist_ok=True)
        allrows = {r['key']: r for r in C.species_rows()}
        for key in sys.argv[2:] or list(rows):
            prot, info = longest_per_gene(allrows[key])
            with open(FAA / f'{key}.faa', 'w') as g:
                for pid, s in prot.items():
                    g.write(f'>{pid}\n{s}\n')
            json.dump(info, open(FAA / f'{key}.json', 'w'))
            print(key, json.dumps(info), flush=True)
    elif cmd == 'run':
        cpus = int(sys.argv[2]) if len(sys.argv) > 2 else 4
        jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 2
        RES.mkdir(parents=True, exist_ok=True)

        def one(key):
            if list((RES / key).glob(f'short_summary*{LINEAGE}*.txt')):
                return key, 0
            p = subprocess.run(PIXI + ['busco', '-i', str(FAA / f'{key}.faa'), '-m', 'prot', '-l', str(FAA / 'busco_downloads' / 'lineages' / LINEAGE), '--offline',
                                '-c', str(cpus), '-o', key, '--out_path', str(RES), '-f'], capture_output=True, text=True)
            (RES / f'{key}.log').write_text(p.stdout[-20000:] + p.stderr[-20000:])
            return key, p.returncode
        with ThreadPoolExecutor(jobs) as ex:
            for key, rc in ex.map(one, list(rows)):
                print(key, 'rc', rc, flush=True)
    elif cmd == 'summary':
        out = {}
        for key in rows:
            js = sorted((RES / key).glob(f'short_summary*{LINEAGE}*.json'))
            if not js:
                out[key] = None
                continue
            r = json.load(open(js[0]))['results']
            out[key] = dict(lineage=LINEAGE, complete=r['Complete percentage'], single=r['Single copy percentage'], duplicated=r['Multi copy percentage'],
                            fragmented=r['Fragmented percentage'], missing=r['Missing percentage'], n=r['n_markers'], one_line=r['one_line_summary'],
                            gate_ok=r['Complete percentage'] >= 97.0)
            print(key, out[key]['one_line'], 'OK' if out[key]['gate_ok'] else 'BELOW 97%')
        json.dump(out, open(RES / 'busco_summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
