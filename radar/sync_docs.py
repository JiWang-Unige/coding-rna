#!/usr/bin/env python3
"""Regenerate the coding-rna dashboard's `docs` collection payload from the project documentation layer.
Usage: python3 radar/sync_docs.py <out_dir>
Writes one JSON per doc (fields project/path/title/size/text/updatedAt) plus batch*.json manifests
(<=50 writes each, <=~400 KB each). Then write them with ArtifactData action "batch" on the dashboard URL
(collection docs, doc_id = JSON basename without .json). When a doc already exists in the collection, pin
its current version with if_version (list the collection first)."""
import json, os, sys, glob, re, datetime
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
files = sorted(glob.glob(f'{ROOT}/docs/review_20261016/[0-9][0-9]_*.md')) + [f'{ROOT}/{p}' for p in [
    'docs/paper_outline_20261021.md', 'docs/paper_outline_20261027.md', 'docs/paper_outline_revision_C1_20261027.md', 'docs/paper_outline_revision_M36_20261101.md', 'docs/REORG_20261016.md', 'docs/README.md', 'AGENTS.md', 'README.md',
    'experiments/M35-PLANTCAD2-LORA/protocol.md', 'experiments/M33-NTV2-FOLLOWUP/protocol.md',
    'radar/seed/SUMMARY.md']]
PER_DOC = 240_000; TOTAL = 1_000_000; BATCH_BYTES = 400_000
w = []; total = 0; batches = [[]]; cur = 0
for f in files:
    if not os.path.exists(f): continue
    rel = os.path.relpath(f, ROOT); t = open(f, encoding='utf-8', errors='replace').read()
    if len(t.encode()) > PER_DOC: t = t[:120_000] + '\n\n[... 截断：原文更长，看本地文件 ...]'
    if total + len(t.encode()) > TOTAL: print('skip (total budget)', rel); continue
    total += len(t.encode())
    did = re.sub(r'[^A-Za-z0-9_\-.]', '_', rel.replace('/', '__'))[:180]
    title = next((l.lstrip('# ').strip()[:120] for l in t.splitlines() if l.strip()), rel)
    doc = dict(project='crna', path=rel, title=title, size=len(t), text=t, updatedAt=datetime.date.today().isoformat())
    p = f'{out}/{did}.json'; json.dump(doc, open(p, 'w'), ensure_ascii=False)
    entry = dict(op='set', collection='docs', doc_id=did, file_path=p)
    sz = os.path.getsize(p)
    if len(batches[-1]) >= 50 or cur + sz > BATCH_BYTES: batches.append([]); cur = 0
    batches[-1].append(entry); cur += sz; w.append(entry)
for i, b in enumerate(batches): json.dump(b, open(f'{out}/batch{i}.json', 'w'), ensure_ascii=False)
print(len(w), 'docs,', total, 'bytes,', len(batches), 'batches')
