"""R1 reference-free run graph and exact ordered interval/path optimization."""
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
import heapq
import itertools
import math
import numpy as np

STOPS = {'TAA', 'TAG', 'TGA'}
STATES = ('',) + tuple('ACGT') + tuple(map(''.join, itertools.product('ACGT', repeat=2)))
AUTOMATON = STATES + ('INIT_A','INIT_AT')  # ATG not completed in a 1/2-bp first CDS.


@dataclass(frozen=True)
class Exon:
    run: int
    a: int
    b: int
    left: str
    right: str
    score: float


@dataclass(frozen=True)
class Path:
    score: float = 0.0
    span: int = 0
    genes: tuple = ()
    active: tuple = ()
    start: int = 0


def genomic_chain(chain, length, strand):
    return tuple(chain) if strand == '+' else tuple((length-b, length-a) for a,b in reversed(chain))


def tie_key(path, length, strand):
    # Compare the eventual completed output after the SAME future continuation.
    # A high sentinel makes an open positive-strand chain's extension sort after
    # an earlier extra exon, instead of treating an unfinished prefix as shorter.
    genes = list(path.genes)
    if path.active:
        genes.append(path.active + (((length+1, length+1),) if strand == '+' else ()))
    return tuple(sorted(genomic_chain(g, length, strand) for g in genes))


def better(a, b, length, strand):
    if a is None:
        return b
    if b is None:
        return a
    if a.score != b.score:
        return a if a.score > b.score else b
    if a.span != b.span:
        return a if a.span < b.span else b
    return a if tie_key(a, length, strand) <= tie_key(b, length, strand) else b


class PrefixBest:
    """Fenwick prefix maxima over original run IDs (strictly increasing)."""
    def __init__(self, n, length, strand):
        self.tree = [None]*(n+1)
        self.length, self.strand = length, strand

    def put(self, index, path):
        i = index+1
        while i < len(self.tree):
            self.tree[i] = better(self.tree[i], path, self.length, self.strand)
            i += i & -i

    def before(self, index):
        result, i = None, index
        while i:
            result = better(result, self.tree[i], self.length, self.strand)
            i -= i & -i
        return result


class ORF:
    def __init__(self, sequence):
        # N/other unknown bases cannot match any character of a stop codon.
        # C also occurs in no stop codon: this preserves stop matching exactly,
        # including exon joins, with the contracted 21 suffix states. Actual
        # sequence and GFF remain unchanged; ATG/terminal motifs use real bases.
        self.sequence = sequence.upper()
        self.s = ''.join(c if c in 'ACGT' else 'C' for c in self.sequence)
        encoded = np.frombuffer(self.s.encode('ascii'), dtype=np.uint8)
        hit = ((encoded[:-2] == 84) & (((encoded[1:-1] == 65) &
               ((encoded[2:] == 65) | (encoded[2:] == 71))) |
               ((encoded[1:-1] == 71) & (encoded[2:] == 65))))
        positions = np.flatnonzero(hit)
        self.stops = [positions[positions % 3 == phase] for phase in range(3)]

    def advance(self, state, exon):
        a,b = exon.a,exon.b
        terminal = exon.right == 'stop'
        initial = state.startswith('INIT_') or exon.left == 'start'
        if state.startswith('INIT_'):
            state = state[5:]
        if exon.left == 'start' and self.sequence[a:a+3] != 'ATG':
            return None
        if terminal and self.sequence[b-3:b] not in STOPS:
            return None
        if terminal and (len(state)+b-a) % 3:
            return None
        # Resolve at most one join-spanning codon, then range-query stops.
        need = (3-len(state)) % 3
        if state and b-a < need:
            suffix = state + self.s[a:b]
            if initial and not 'ATG'.startswith(suffix):
                return None
            return None if terminal else ('INIT_'+suffix if initial else suffix)
        if state:
            codon = state+self.s[a:a+need]
            if initial and codon != 'ATG':
                return None
            if codon in STOPS:
                return '' if terminal and a+need == b else None
            a += need
            initial = False
            if terminal and a == b:
                return None
        if initial:
            if b-a < 3:
                return None if terminal else 'INIT_'+self.s[a:b]
            if self.s[a:a+3] != 'ATG':
                return None
        end = b-3 if terminal else b
        if end < a:
            return None
        positions = self.stops[a % 3]
        i = bisect_left(positions, a)
        if i < len(positions) and positions[i]+3 <= end:
            return None
        trailing = (end-a) % 3
        suffix = self.s[end-trailing:end] if trailing else ''
        return None if terminal and suffix else suffix


def solve(sequence, exons, intron_prefix, strand='+'):
    """Coordinate sweep; finished genes and splice donors activate separately.

    Intron edges require both end+4<=acceptor and previous_run<current_run.
    Prefix maxima handle run order even when +/-6 motifs reorder coordinates.
    Donor scores subtract G-prefix; span subtracts the open gene's start.
    """
    n = len(sequence)
    exons = sorted(exons, key=lambda e: (e.a,e.b,e.run,e.left,e.right))
    if not exons:
        return Path()
    orf = ORF(sequence)
    donors = {s: PrefixBest(max(e.run for e in exons)+1,n,strand) for s in AUTOMATON}
    pending, serial = [], itertools.count()
    complete = Path()
    for exon in exons:
        while pending and pending[0][0] <= exon.a:
            _at, _id, kind, state, run, p = heapq.heappop(pending)
            if kind == 'gene':
                complete = better(complete,p,n,strand)
            else:
                donors[state].put(run,p)
        if exon.left == 'start':
            incoming = [('',Path(complete.score,complete.span-exon.a,complete.genes,(),exon.a))]
        else:
            incoming = [(s, donors[s].before(exon.run)) for s in AUTOMATON]
        for state,p in incoming:
            if p is None:
                continue
            following = orf.advance(state,exon)
            if following is None:
                continue
            value = p.score + exon.score
            if exon.left == 'acceptor':
                value += float(intron_prefix[exon.a])
            chain = p.active + ((exon.a,exon.b),)
            if exon.right == 'stop':
                if sum(b-a for a,b in chain) < 6:
                    continue
                candidate = Path(value,p.span+exon.b,p.genes+(chain,))
                heapq.heappush(pending,(exon.b,next(serial),'gene','',exon.run,candidate))
            else:
                candidate = Path(value-float(intron_prefix[exon.b]),p.span,p.genes,chain,p.start)
                heapq.heappush(pending,(exon.b+4,next(serial),'donor',following,exon.run,candidate))
    for _at,_id,kind,_state,_run,p in pending:
        if kind == 'gene':
            complete = better(complete,p,n,strand)
    return complete


def prefixes(region):
    return tuple(np.r_[0., np.cumsum(region[:,c].astype(np.float64)-region[:,0].astype(np.float64))]
                 for c in (1,2))


def candidates(sequence, region, boundary):
    # No reference, saved trace, oracle manifest, or species identity is read.
    import redecode_error_decomposition as D
    states = D.m25.region_state_path(region,0.4)
    cds_prefix,intron_prefix = prefixes(region)
    changes = np.r_[0,np.flatnonzero(states[1:] != states[:-1])+1,len(states)]
    runs = [(int(a),int(b)) for a,b in zip(changes,changes[1:]) if states[a] == D.m25.C]
    exons=[]
    thresholds={'start':0.5,'stop':0.5,'donor':0.1,'acceptor':0.1}
    for rid,(a,b) in enumerate(runs):
        options={name:set() for name in thresholds}
        for position in (a,b):
            if not 0 < position < len(states):
                continue
            for event in D.m25.TRANSITIONS.get((int(states[position-1]),int(states[position])),()):
                anchor=D.event_anchor(position,event)
                channel=D.m25.BOUNDARY_NAMES.index(event)
                cutoff=math.log(thresholds[event]/(1-thresholds[event]))
                for pos in D.motif_positions(sequence,event,anchor):
                    if float(boundary[pos,channel]) >= cutoff:
                        options[event].add(pos)
        for left,right in itertools.product(('start','acceptor'),('donor','stop')):
            li,ri=(D.m25.BOUNDARY_NAMES.index(x) for x in (left,right))
            for p,q in itertools.product(sorted(options[left]),sorted(options[right])):
                start=p+(2 if left=='acceptor' else 0)
                end=q+(3 if right=='stop' else 0)
                if start >= end:
                    continue
                score=float(cds_prefix[end]-cds_prefix[start])
                for pos,channel,event in ((p,li,left),(q,ri,right)):
                    score += float(boundary[pos,channel])-math.log(thresholds[event]/(1-thresholds[event]))
                exons.append(Exon(rid,start,end,left,right,score))
    return exons,intron_prefix,{'original_CDS_runs':len(runs),'candidate_CDS_fragments':len(exons)}


def models(path):
    result=[]
    for chain in path.genes:
        phases,total=[],0
        for a,b in chain:
            phases.append((3-total%3)%3)
            total += b-a
        result.append({'cds':list(chain),'phase':phases,
                       'start_codon':(chain[0][0],chain[0][0]+3),
                       'stop_codon':(chain[-1][1]-3,chain[-1][1]),
                       'region_span':(chain[0][0],chain[-1][1]),
                       'boundary_scores':{}})
    return result
