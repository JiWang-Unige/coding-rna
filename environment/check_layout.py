#!/usr/bin/env python3
"""Static layout checks for the reorganised repository (stdlib only, no torch needed).

Run from the repository root:  python3 environment/check_layout.py

1. every tracked/untracked .py compiles (syntax check without writing .pyc)
2. every import in pipeline/m28 and tests/ resolves (third-party and stdlib excluded)
3. every script / test referenced by pipeline/m28/sbatch exists and uses the new layout
4. every relative Markdown/HTML link resolves (links into ignored runtime dirs such as outputs/ are only counted)
Exit status is non-zero when anything is broken.
"""
import ast
import importlib.util
import posixpath
import re
import subprocess
import sys
import sysconfig
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
PIPE = ROOT / "pipeline" / "m28"
THIRD = {
    "numpy", "torch", "scipy", "pandas", "yaml", "transformers", "tokenizers", "peft", "accelerate", "numba",
    "pyfaidx", "Bio", "einops", "safetensors", "huggingface_hub", "pytest", "tiberius", "jax", "haiku",
    "nucleotide_transformer", "matplotlib", "sklearn", "tqdm", "h5py", "psutil",
}
RUNTIME_TOP = ("outputs", "data", "logs", "refs", "runs")
failures = []


def tracked_files():
    out = subprocess.check_output(["git", "ls-files", "-co", "--exclude-standard"], cwd=ROOT, text=True)
    return out.splitlines()


def is_stdlib(top):
    if top in sys.builtin_module_names:
        return True
    try:
        spec = importlib.util.find_spec(top)
    except Exception:
        return False
    return bool(spec and spec.origin and spec.origin.startswith(sysconfig.get_paths()["stdlib"])
                and "site-packages" not in spec.origin)


def local_exists(name, roots):
    parts = name.split(".")
    return any((r.joinpath(*parts).with_suffix(".py").exists() or r.joinpath(*parts).is_dir()) for r in roots)


def check_compile(files):
    bad = []
    for f in files:
        if f.endswith(".py"):
            try:
                compile((ROOT / f).read_text(encoding="utf-8"), f, "exec")
            except SyntaxError as exc:
                bad.append((f, str(exc)))
    n = sum(f.endswith(".py") for f in files)
    print(f"[1] syntax: {n} python files, {len(bad)} failed")
    failures.extend(bad)


def check_imports():
    roots = [PIPE, PIPE / "scripts", ROOT / "tests"]
    unresolved, total = [], 0
    for base in (PIPE, ROOT / "tests"):
        for path in sorted(base.rglob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                names = []
                if isinstance(node, ast.Import):
                    names = [a.name for a in node.names]
                elif isinstance(node, ast.ImportFrom):
                    if node.level:
                        total += 1
                        target = path.parent
                        for _ in range(node.level - 1):
                            target = target.parent
                        mod = (node.module or "").replace(".", "/")
                        if mod and not ((target / (mod + ".py")).exists() or (target / mod).is_dir()):
                            unresolved.append((str(path.relative_to(ROOT)), node.lineno, "." * node.level + (node.module or "")))
                        continue
                    names = [node.module]
                for name in names:
                    total += 1
                    top = name.split(".")[0]
                    if top in THIRD or is_stdlib(top) or local_exists(name, roots):
                        continue
                    unresolved.append((str(path.relative_to(ROOT)), node.lineno, name))
    print(f"[2] imports: {total} import statements in pipeline/m28 + tests, {len(unresolved)} unresolved")
    failures.extend(unresolved)


def check_sbatch():
    bad, refs, count = [], 0, 0
    for f in sorted((PIPE / "sbatch").rglob("*.sbatch")):
        count += 1
        text = f.read_text()
        rel = str(f.relative_to(ROOT))
        for m in re.findall(r"pipeline/m28/scripts/[A-Za-z0-9_./-]+\.py|tests/[A-Za-z0-9_./-]+\.py", text):
            refs += 1
            if not (ROOT / m).is_file():
                bad.append((rel, m))
        for m in re.findall(r"-p\s+(test_[A-Za-z0-9_]+\.py)|pattern=\s*'(test_[A-Za-z0-9_]+\.py)'", text):
            name = m[0] or m[1]
            refs += 1
            if not (ROOT / "tests" / name).is_file():
                bad.append((rel, "tests/" + name))
        for m in re.findall(r"src\.m28\.([A-Za-z0-9_]+)", text):
            refs += 1
            if not (PIPE / "src" / "m28" / (m + ".py")).is_file():
                bad.append((rel, "src.m28." + m))
        if "scripts/experiments" in text or re.search(r'PYTHONPATH="\$PWD(:|")', text):
            bad.append((rel, "old layout (scripts/experiments or PYTHONPATH=$PWD)"))
    print(f"[3] sbatch: {count} pipeline sbatch files, {refs} path references, {len(bad)} bad")
    failures.extend(bad)


LINK = re.compile(r'(?:!?\[[^\]]*\]\(<?([^)\s>]+)>?(?:\s+"[^"]*")?\))|(?:^\s{0,3}\[[^\]]+\]:\s*(\S+))|(?:(?:href|src)="([^"]+)")', re.M)


def check_links(files):
    fs, dirs = set(files), set()
    for f in files:
        d = posixpath.dirname(f)
        while d:
            dirs.add(d)
            d = posixpath.dirname(d)
    total = ok = runtime = 0
    bad = []
    for f in files:
        if not f.endswith((".md", ".html")) or f.startswith("archive/backups/"):
            continue
        text = (ROOT / f).read_text(encoding="utf-8", errors="replace")
        for m in LINK.finditer(text):
            target = next(g for g in m.groups() if g)
            if re.match(r"^([a-zA-Z][a-zA-Z0-9+.\-]*:|#|/)", target):
                continue
            path = unquote(target.split("#", 1)[0])
            if not path:
                continue
            total += 1
            res = posixpath.normpath(posixpath.join(posixpath.dirname(f), path))
            if res in fs or res in dirs or res == ".":
                ok += 1
            elif res.split("/")[0] in RUNTIME_TOP:
                runtime += 1
            else:
                bad.append((f, target))
    print(f"[4] links: {total} relative links, {ok} resolve, {runtime} point into ignored runtime dirs, {len(bad)} broken")
    failures.extend(bad)


if __name__ == "__main__":
    files = tracked_files()
    check_compile(files)
    check_imports()
    check_sbatch()
    check_links(files)
    for item in failures:
        print("FAIL", item)
    sys.exit(1 if failures else 0)
