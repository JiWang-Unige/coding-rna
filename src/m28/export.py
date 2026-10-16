"""Primary CDS-chain GFF3 export; no UTR or alternative-isoform claim."""
def as_transcript(chain,seqid,strand,identifier,score=None):
    parts=[tuple(p) for p in chain]
    phases={};offset=0
    for a,b in parts if strand=='+' else parts[::-1]:
        phases[a,b]=str((3-offset%3)%3)
        offset+=b-a
    return {'id':identifier+'.t1','gene_id':identifier,'seqid':seqid,'strand':strand,
            'CDS':[(a,b,phases[a,b]) for a,b in parts],'exon':parts,'partial':False,'score':score}

def write_gff3(path,transcripts,lengths,source):
    with path.open('x') as f:
        f.write('##gff-version 3\n')
        for seqid,length in lengths.items(): f.write(f'##sequence-region {seqid} 1 {length}\n')
        for tx in transcripts:
            seqid,strand=tx['seqid'],tx['strand']
            a,b=tx['CDS'][0][0],tx['CDS'][-1][1]
            score='.' if tx.get('score') is None else format(tx['score'],'.8g')
            gene,tid=tx['gene_id'],tx['id']
            def row(kind,x,y,phase,attrs):
                f.write('\t'.join(map(str,(seqid,source,kind,x+1,y,score,strand,phase,attrs)))+'\n')
            row('gene',a,b,'.',f'ID={gene};gene_biotype=protein_coding')
            row('mRNA',a,b,'.',f'ID={tid};Parent={gene}')
            for i,(x,y,phase) in enumerate(tx['CDS'],1):
                row('exon',x,y,'.',f'ID={tid}.exon{i};Parent={tid}')
                row('CDS',x,y,phase,f'ID={tid}.cds{i};Parent={tid}')
