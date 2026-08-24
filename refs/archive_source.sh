#!/usr/bin/env bash
# archive_source.sh — best-effort archive of a paper + its code repo into refs/.
#
# Used by /sota-inventory (--type sota) and /note-add (--type note).
# Best-effort: download failures are recorded, never fatal — so the caller can
# still proceed and the user can drop the PDF in manually later.
#
# Usage:
#   bash refs/archive_source.sh --slug <slug> [--arxiv <id>] [--pdf-url <url>] \
#        [--repo <git-url>] [--title "..."] [--type sota|note] [--why "..."] [--refs-dir refs]
set -uo pipefail

SLUG="" ARXIV="" PDF_URL="" REPO="" TITLE="" TYPE="note" WHY="" REFS="refs"
while [ $# -gt 0 ]; do
  case "$1" in
    --slug) SLUG="$2"; shift 2;;
    --arxiv) ARXIV="$2"; shift 2;;
    --pdf-url) PDF_URL="$2"; shift 2;;
    --repo) REPO="$2"; shift 2;;
    --title) TITLE="$2"; shift 2;;
    --type) TYPE="$2"; shift 2;;
    --why) WHY="$2"; shift 2;;
    --refs-dir) REFS="$2"; shift 2;;
    *) echo "WARN: unknown arg $1" >&2; shift;;
  esac
done

if [ -z "$SLUG" ]; then echo "ERROR: --slug required" >&2; exit 2; fi
mkdir -p "$REFS/pdfs" "$REFS/repos" "$REFS/dossiers"

PDF_STATUS="none"; REPO_STATUS="none"; REPO_COMMIT=""

# --- 1. PDF (arXiv id preferred, else explicit url) ---
PDF_PATH="$REFS/pdfs/${SLUG}.pdf"
DL_URL=""
[ -n "$ARXIV" ] && DL_URL="https://arxiv.org/pdf/${ARXIV}.pdf"
[ -z "$DL_URL" ] && [ -n "$PDF_URL" ] && DL_URL="$PDF_URL"
if [ -n "$DL_URL" ]; then
  if [ -s "$PDF_PATH" ]; then
    PDF_STATUS="exists"
  elif curl -fsSL --max-time 60 "$DL_URL" -o "$PDF_PATH" 2>/dev/null && [ -s "$PDF_PATH" ]; then
    PDF_STATUS="downloaded"
  else
    rm -f "$PDF_PATH"
    PDF_STATUS="failed($DL_URL)"
  fi
fi

# --- 2. Repo (shallow clone; record commit) ---
REPO_DIR="$REFS/repos/${SLUG}"
if [ -n "$REPO" ]; then
  if [ -d "$REPO_DIR/.git" ]; then
    REPO_STATUS="exists"
    REPO_COMMIT="$(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null || echo '?')"
  elif git clone --depth 1 --quiet "$REPO" "$REPO_DIR" 2>/dev/null; then
    REPO_STATUS="cloned"
    REPO_COMMIT="$(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null || echo '?')"
  else
    # too big / private / no git access — keep a link stub instead
    printf '# %s (repo not cloned)\n\n- URL: %s\n- Reason: clone failed (private/large/no-net). Clone manually if needed.\n' \
      "$SLUG" "$REPO" > "$REFS/repos/${SLUG}.link.md"
    REPO_STATUS="link-only"
  fi
fi

# --- 3. Dossier skeleton (only if missing — never clobber filled-in detail) ---
DOSSIER="$REFS/dossiers/${SLUG}.md"
if [ ! -f "$DOSSIER" ]; then
  {
    echo "# Dossier: ${TITLE:-$SLUG}"
    echo
    echo "- slug: \`$SLUG\` · type: $TYPE · added: $(date +%F)"
    echo "- Links: ${ARXIV:+arXiv:$ARXIV } ${REPO:+repo:$REPO}"
    echo "- PDF: refs/pdfs/${SLUG}.pdf ($PDF_STATUS)"
    echo "- Repo: refs/repos/${SLUG}/ ($REPO_STATUS${REPO_COMMIT:+ @ $REPO_COMMIT})"
    [ -n "$WHY" ] && echo "- Why relevant: $WHY"
    echo
    echo "## Dataset source (⏳ verify via WebFetch)"
    echo "- Which dataset / version / where obtained:"
    echo "- Public download / license:"
    echo
    echo "## Metric implementation (⏳ verify)"
    echo "- Metric name + exact definition (e.g., segment F1 boundary rule):"
    echo "- Official impl / script / library:"
    echo
    echo "## Split scheme (⏳ verify)"
    echo "- Train/val/test split source + leakage notes:"
    echo
    echo "## Weights / license"
    echo "- Pretrained weights URL + version + license:"
    echo
    echo "## Reproducibility notes"
    echo "- Setup difficulty / known gotchas:"
    echo
    echo "## Relevance to our project"
    echo "- ${WHY:-}"
  } > "$DOSSIER"
  DOSSIER_STATUS="created"
else
  DOSSIER_STATUS="exists"
fi

# --- 4. Append to index (avoid dup slug row) ---
INDEX="$REFS/sources.md"
[ -f "$INDEX" ] || printf '# Archived Sources Index\n\n| slug | title | type | pdf | repo | dossier | added_by | date |\n|---|---|---|---|---|---|---|---|\n' > "$INDEX"
if ! grep -q "| \`\?${SLUG}\`\? |" "$INDEX" 2>/dev/null && ! grep -q "^| ${SLUG} |" "$INDEX" 2>/dev/null; then
  printf '| %s | %s | %s | %s | %s | %s | %s | %s |\n' \
    "$SLUG" "${TITLE:-}" "$TYPE" "$PDF_STATUS" "$REPO_STATUS" "dossiers/${SLUG}.md" "${ADDED_BY:-archive_source}" "$(date +%F)" >> "$INDEX"
fi

echo "ARCHIVED slug=$SLUG pdf=$PDF_STATUS repo=$REPO_STATUS${REPO_COMMIT:+@$REPO_COMMIT} dossier=$DOSSIER_STATUS"
