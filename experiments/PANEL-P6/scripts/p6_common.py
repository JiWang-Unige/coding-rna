"""PANEL-P6 shared setup (stdlib only): the P6 panel baselines of paper_outline_20261027 P6 on the core-10 independent panel (report 25 §9), scored against our P5 T2 run.
Protocol / record: docs/review_20261016/31_P6_panel_baselines.md.  C1 criterion: report 28 §6; seed rule: report 28 §7.

Paths on ont: code ~/coding-rna-reorg/experiments/PANEL-P6, outputs ~/coding-rna/outputs/PANEL-P6, panel data ~/coding-rna/data/panel (PANEL-PREP, frozen).
"""
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
HOME = Path(os.path.expanduser('~'))
DATA = HOME / 'coding-rna' / 'data' / 'panel'
RAW = DATA / 'raw'
PREP = DATA / 'prepared'
OUT = Path(os.environ.get('P6_OUT', HOME / 'coding-rna' / 'outputs' / 'PANEL-P6'))
GENOMES = OUT / 'genomes'
T2_PANEL = HOME / 'coding-rna' / 'outputs' / 'M34-NTV2-SCALEUP' / 'panel' / 'T2'      # P5: <key>/<contig>/prediction/full_grammar.gff3
NATIVE = OUT / 'native'          # native/<Helixer|Tiberius|ANNEVO>/<key>/{prediction.gff3, DONE}
CARBON = OUT / 'carbonA'         # carbonA/<key>/<input stem>/{*.raw_probs.parquet, prediction.gff}, carbonA/<key>/DONE
VIP = OUT / 'vipsania'           # vipsania/runs/<key>/<ft|pre>/{vipsania.gff3, <mode>_DONE}
RUNS = OUT / 'runs'              # converted predictions: runs/<tool>/<key>/<contig>/prediction/full_grammar.gff3 (m34_infer layout)
RESULTS = OUT / 'results'

ALL12 = ['dryas_octopetala', 'eschscholzia_californica', 'akebia_trifoliata', 'hakea_prostrata', 'chenopodium_album', 'zizania_latifolia', 'themeda_triandra',
         'juncus_effusus', 'rhynchospora_breviuscula', 'acorus_gramineus', 'dipteris_shenzhenensis', 'amborella_trichopoda']
CORE10 = [k for k in ALL12 if k not in ('chenopodium_album', 'dipteris_shenzhenensis')]       # same list and order as m34_panel_common.CORE10


def species():
    rows = json.load(open(EXP / 'PANEL-PREP' / 'scripts' / 'panel_species.json'))
    return {r['key']: r for r in rows}


def raw_dir(key):
    r = species()[key]
    return RAW / f"{r['accession']}_{r['assembly_name']}"


def lengths(key):
    """main (scored) sequences of a panel species: {contig: bp} from PANEL-PREP prepared/<key>/lengths.json"""
    return {k: int(v) for k, v in json.load(open(PREP / key / 'lengths.json')).items()}


def by_size(keys):
    """smallest whole assembly first (NCBI total, panel_species.json doc10_total_mb)"""
    sp = species()
    return sorted(keys, key=lambda k: (sp[k]['doc10_total_mb'], k))
