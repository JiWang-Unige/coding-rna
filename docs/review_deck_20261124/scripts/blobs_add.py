"""Record uploaded figure assets: python blobs_add.py <id> <figure png> [<id> <png> ...]
Writes blobs.json ({id: {path, sha256}}) next to this folder's figures/ and a flat copy
blobs_flat.json ({id: path}) for render_preview.py, which only reads the flat form."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
p = os.path.join(HERE, "blobs.json")
b = json.load(open(p)) if os.path.exists(p) else {}
args = sys.argv[1:]
for i, f in zip(args[::2], args[1::2]):
    f = os.path.abspath(f)
    b[i] = {"path": f, "sha256": hashlib.sha256(open(f, "rb").read()).hexdigest()}
json.dump(b, open(p, "w"), indent=1)
json.dump({k: v["path"] for k, v in b.items()}, open(os.path.join(HERE, "blobs_flat.json"), "w"), indent=1)
print(len(b), "assets recorded")
