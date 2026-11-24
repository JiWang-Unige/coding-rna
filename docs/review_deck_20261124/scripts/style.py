"""Shared figure style and provenance checks for the review decks (2 Nov 2026, updated 24 Nov 2026).

Slides are 1920x1080 px and figures are placed at 100 px per inch, so a figure
made with figure(1664, 560) fills a 1664x560 px slot exactly. Saved at dpi=200
(2x pixel density). Tick labels land at ~17-18 px on the slide.

cite(path, line, *tokens) is the provenance check every figure uses: it asserts
that each token occurs verbatim on that line of the repository file (the result
report or table the number comes from) and returns the tokens as floats. When a
report is edited and a number moves or changes, the figure script stops.
"""
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402

for _f in (
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
):
    if os.path.exists(_f):
        fm.fontManager.addfont(_f)

REPO = "/Users/jiwang/Desktop/coding-rna-reorg"
OUTDIR = os.path.join(REPO, "docs/review_deck_20261124/figures")

# Slot sizes on the slide (px)
SLOT_WIDE = (1664, 560)
SLOT_LEFT = (1060, 560)

# Palette (matches the deck)
BG = "#FAF9F5"           # slide background
INK = "#1B2430"          # headings, axis labels
BODY = "#3D4757"         # tick labels, annotations
MUTED = "#6B7380"        # spines, secondary text
GRID = "#E4E1D8"         # light reference lines
BLUE = "#1F5F8B"         # our final recipe (PlantCAD2 + LoRA)
BLUE_LIGHT = "#9DBBD6"   # our other variants
ORANGE = "#C2662D"       # Carbon-A, or the contrast arm (fill-background)
ORANGE_LIGHT = "#E8B48F"
GRAY = "#A3A9B2"         # other tools, references
GRAY_LIGHT = "#D3D6DA"
SLATE = "#5B6B7D"        # GENERanno (the 16 October model)


def apply():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 18,
        "axes.titlesize": 19,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.titlepad": 14,
        "axes.labelsize": 18,
        "xtick.labelsize": 17,
        "ytick.labelsize": 17,
        "legend.fontsize": 16,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "xtick.color": BODY,
        "ytick.color": BODY,
        "text.color": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "axes.facecolor": BG,
        "figure.facecolor": BG,
        "savefig.facecolor": BG,
        "axes.grid": False,
        "legend.frameon": False,
        "lines.linewidth": 2.4,
        "lines.markersize": 9,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })


def figure(width_px, height_px, **kwargs):
    """Figure that fills a width_px x height_px slot on the slide."""
    apply()
    return plt.subplots(figsize=(width_px / 100, height_px / 100), **kwargs)


def axes_px(fig, left, top, width, height):
    """Axes placed in slot pixels (origin top-left)."""
    W, H = fig.get_size_inches() * 100
    return fig.add_axes([left / W, 1 - (top + height) / H, width / W, height / H])


def legend_row(fig, left, top, items, gap=36, fontsize=17):
    """Hand-drawn legend in one row at slot px (left, top). items: (kind, face, edge, text);
    kind is a matplotlib marker ('o', 'D', 's') or 'line' / 'dash'. Widths are measured."""
    W, H = fig.get_size_inches() * 100
    ax = axes_px(fig, left, top, W - left, 40)
    ax.axis("off")
    ax.set_xlim(0, W - left)
    ax.set_ylim(0, 40)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    x = 10
    for kind, fc, ec, text in items:
        if kind in ("line", "dash"):
            ax.plot([x - 10, x + 14], [20, 20], color=ec, lw=3, ls="--" if kind == "dash" else "-",
                    clip_on=False)
            x += 8
        else:
            ax.scatter([x], [20], marker=kind, s=110 if kind == "o" else 70, facecolor=fc,
                       edgecolor=ec, linewidth=2.6, clip_on=False)
        t = ax.text(x + 18, 20, text, va="center", fontsize=fontsize, color=BODY)
        wpx = t.get_window_extent(renderer=r).width / fig.dpi * 100
        x += 18 + wpx + gap
    return ax


def save(fig, name):
    """Write <name>.png (dpi 200) to OUTDIR; keeps the exact slot size."""
    os.makedirs(OUTDIR, exist_ok=True)
    png = os.path.join(OUTDIR, name + ".png")
    fig.savefig(png, dpi=200)
    plt.close(fig)
    return png


def _num(tok):
    t = tok.replace("−", "-").replace("+", "").replace(",", "").replace("%", "")
    return float(re.sub(r"[^0-9.\-]", "", t))


def line(path, lineno):
    with open(os.path.join(REPO, path), encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    return lines[lineno - 1]


def cite(path, lineno, *tokens):
    """Assert each token occurs on path:lineno; return the tokens as floats."""
    text = line(path, lineno)
    for t in tokens:
        assert t in text, f"{path}:{lineno} does not contain {t!r}\n  line: {text[:240]}"
    return [_num(t) for t in tokens]


def has(path, lineno, *phrases):
    """Assert each phrase occurs on path:lineno (no number conversion)."""
    text = line(path, lineno)
    for t in phrases:
        assert t in text, f"{path}:{lineno} does not contain {t!r}\n  line: {text[:240]}"


def cells(path, lineno):
    """Cells of a Markdown table row, with ** removed."""
    text = line(path, lineno).replace("**", "")
    assert text.startswith("|"), f"{path}:{lineno} is not a table row"
    return [c.strip() for c in text.strip().strip("|").split("|")]


def close(a, b, tol=0.006):
    assert abs(a - b) <= tol, f"{a} != {b}"


def find(path, *phrases):
    """Line number (1-based) of the single line of path that contains every phrase."""
    with open(os.path.join(REPO, path), encoding="utf-8") as fh:
        hits = [i + 1 for i, l in enumerate(fh.read().splitlines()) if all(p in l for p in phrases)]
    assert len(hits) == 1, f"{path}: {phrases} found on lines {hits}"
    return hits[0]
