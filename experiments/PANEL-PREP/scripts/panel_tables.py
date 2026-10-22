"""Aggregate the ont-side products into the experiment's results tables (run on ont after prep/stats/homology; copy results/ back to the Mac):
  download_manifest.tsv  (accession, files, md5 local vs NCBI, download time, NCBI annotation name/date/method, BUSCO)
  species_table.tsv      (tier, kept sequences / bp, dropped, genes, CDS-assessable references, exon strata, tool exposure)
  reference_stats.tsv    (full per-species reference statistics)
  disk_usage.txt
usage: panel_tables.py <results_dir>
"""
import json
import subprocess
import sys
from pathlib import Path

import panel_common as P


def main(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    tools = json.load(open(P.HERE / 'tool_exposure.json'))
    man, spt, rst = [], [], []
    for key, row in P.species().items():
        raw = P.raw_dir(row)
        md5 = {}
        for line in open(raw / 'local_md5.txt'):
            f = line.split()
            md5[f[1]] = (f[0], f[2].split('=')[1], f[-1])
        dl = open(raw / 'download_date.txt').readline().strip()
        sp = json.load(open(P.PREP / key / 'splits.json')) if (P.PREP / key / 'splits.json').exists() else {}
        nc = sp.get('ncbi') or {}
        st = json.load(open(P.PREP / key / 'reference_stats.json')) if (P.PREP / key / 'reference_stats.json').exists() else {}
        sizes = {s: (raw / (P.stem(row) + s)).stat().st_size / 1e6 for s in ('_genomic.fna.gz', '_genomic.gff.gz', '_protein.faa.gz')}
        for s in ('_genomic.fna.gz', '_genomic.gff.gz', '_protein.faa.gz', '_assembly_report.txt', '_assembly_stats.txt'):
            f = P.stem(row) + s
            man.append([key, row['accession'], row['assembly_name'], f, f'{(raw / f).stat().st_size / 1e6:.1f}', md5[f][0], md5[f][1], md5[f][2], dl,
                        nc.get('annotation_name'), nc.get('annotation_release_date'), nc.get('annotation_provider'), nc.get('annotation_method'),
                        (nc.get('busco') or {}).get('lineage'), (nc.get('busco') or {}).get('complete')])
        gc = nc.get('gene_counts') or {}
        spt.append([key, row['organism'], row['tier'], row['source'], row['accession'], row['assembly_name'], nc.get('annotation_name'), nc.get('annotation_release_date'),
                    nc.get('assembly_level'), f"{sizes['_genomic.fna.gz']:.0f}/{sizes['_genomic.gff.gz']:.1f}/{sizes['_protein.faa.gz']:.1f}",
                    nc.get('total_sequence_length'), sp.get('rule', '').split(' (')[0], len(sp.get('test', [])), sp.get('test_bp'), sp.get('dropped_sequences'), sp.get('dropped_bp'),
                    sp.get('dropped_annotated_genes'), gc.get('total'), gc.get('protein_coding'), sp.get('gff_genes'), (sp.get('gene_biotypes') or {}).get('protein_coding'),
                    st.get('coding_genes_parsed'), st.get('cds_assessable_references'), st.get('complete_supported_nonoverlapping'),
                    '/'.join(str(st['exon_hist'][k]) for k in ('1', '2-5', '>=6')) if st else None, st.get('median_exons'), st.get('spans_longer_than_window'),
                    f"{st['nonstandard_splice_fraction']:.4f}" if st.get('nonstandard_splice_fraction') is not None else None,
                    (nc.get('busco') or {}).get('complete'), md5[P.stem(row) + '_genomic.gff.gz'][0][:10],
                    tools[key]['carbon_a'], tools[key]['helixer_a0080'], tools[key]['tiberius_angiosperms'], tools[key]['annevo_paper'], tools[key]['plantgeneann'], tools[key]['plantcad2']])
        if st:
            rst.append([key, st['contigs'], st['bp'], st['coding_genes_parsed'], st['primary_coding_genes'], st['cds_assessable_references'], f"{st['assessable_per_Mb']:.1f}",
                        st['complete_supported_nonoverlapping'], st['exon_hist']['1'], st['exon_hist']['2-5'], st['exon_hist']['>=6'], st['median_exons'], st['median_cds_bp'],
                        json.dumps(st['cds_bp_hist']), st['spans_longer_than_window'], f"{st['nonstandard_splice_fraction']:.4f}", json.dumps(st['unsupported_reasons'])])
    heads = dict(download_manifest=['key', 'accession', 'assembly_name', 'file', 'MB', 'md5_local', 'md5_ncbi', 'status', 'download_time', 'ncbi_annotation_name', 'ncbi_annotation_date',
                                    'annotation_provider', 'annotation_method', 'busco_lineage', 'busco_complete'],
                 species_table=['key', 'organism', 'tier', 'source', 'accession', 'assembly_name', 'annotation_name', 'annotation_date', 'assembly_level', 'gz_MB_fna/gff/faa',
                                'ncbi_total_bp', 'sequence_rule', 'kept_seqs', 'kept_bp', 'dropped_seqs', 'dropped_bp', 'dropped_genes', 'ncbi_genes_total', 'ncbi_protein_coding',
                                'kept_gff_genes', 'kept_protein_coding_genes', 'parsed_coding_genes', 'cds_assessable_refs', 'complete_supported_nonoverlapping', 'exons_1/2-5/>=6',
                                'median_exons', 'cds_span_gt_24576', 'nonstd_splice_frac', 'ref_busco_complete', 'gff_md5_10', 'exp_carbon_a', 'exp_helixer', 'exp_tiberius', 'exp_annevo', 'exp_plantgeneann', 'exp_plantcad2'],
                 reference_stats=['key', 'contigs', 'bp', 'coding_genes_parsed', 'primary_coding_genes', 'cds_assessable_refs', 'assessable_per_Mb', 'complete_supported_nonoverlapping',
                                  'exons_1', 'exons_2-5', 'exons_ge6', 'median_exons', 'median_cds_bp', 'cds_bp_hist', 'cds_span_gt_24576', 'nonstd_splice_frac', 'unsupported_reasons'])
    for name, rows in (('download_manifest', man), ('species_table', spt), ('reference_stats', rst)):
        with open(out / f'{name}.tsv', 'w') as f:
            f.write('\t'.join(heads[name]) + '\n')
            for r in rows:
                f.write('\t'.join('' if x is None else str(x) for x in r) + '\n')
    du = subprocess.run(['bash', '-c', f'du -sh {P.DATA}/* {P.DATA}; df -h {P.DATA} | tail -1'], capture_output=True, text=True).stdout
    (out / 'disk_usage.txt').write_text(du)
    print(du)
    print('rows', len(man), len(spt), len(rst))


if __name__ == '__main__':
    main(sys.argv[1])
