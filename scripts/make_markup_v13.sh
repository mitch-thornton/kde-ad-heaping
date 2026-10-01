#!/usr/bin/env bash
# make_markup_v13.sh
#
# Builds the marked-up manuscript, submission artifact E-3: the version the reviewers read
# last round against the current revision, with every change shown inline.
#
# You use zsh and this is a bash script, so invoke it as
#
#     bash scripts/make_markup_v13.sh
#
# Needs latexdiff on PATH. It ships with TeX Live, so `tlmgr install latexdiff` if the
# command is missing, or `brew install latexdiff`.
#
# Two preparations matter and neither is optional.
#
# First, the current source reads its numbers from tables/numbers_v10_7.tex and its tables
# from tables/*.tex. Diffing it raw would show a reviewer "\figIseBreakRatio" where the old
# version showed a number. So the macros are expanded and the table files inlined before
# the comparison, and the markup then shows real values on both sides.
#
# Second, the wlscirep class declares the abstract before \begin{document}, so latexdiff
# treats it as preamble and leaves it unmarked. The abstract carries the two corrections
# the editor asked for, C1 and C2, so it is diffed separately and spliced back in.
#
# Environment:
#   OLD       the reviewed .tex              (default: ../reviewed/heaped_kde_sr_v12_5.tex)
#   OLDTABLES the table set that shipped with OLD (default: ../reviewed/tables_v12_5)
#   NEW       the current .tex                (default: ../heaped_kde_sr_v13_1.tex)
#   OUTDIR    where the markup is written     (default: ../markup)
#
# This round changes two sentences in Methods, one word in the supplementary-table range,
# and nothing else. A markup that shows more than that is wrong.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
## The round-five markup is against the version the reviewers read, v12.5, not against the
## version that entered review. Set OLD to reviewed/heaped_kde_sr_v9_4.tex, with OLDTABLES
## unset, for the cumulative comparison instead.
OLD="${OLD:-$ROOT/reviewed/heaped_kde_sr_v12_5.tex}"
NEW="${NEW:-$ROOT/heaped_kde_sr_v13_1.tex}"
OUTDIR="${OUTDIR:-$ROOT/markup}"
STEM="heaped_kde_sr_markup"

die() { echo "error: $*" >&2; exit 1; }
command -v latexdiff >/dev/null 2>&1 || die "latexdiff is not on PATH (tlmgr install latexdiff)"
command -v pdflatex  >/dev/null 2>&1 || die "pdflatex is not on PATH"
[ -f "$OLD" ] || die "reviewed source not found at $OLD"
[ -f "$NEW" ] || die "current source not found at $NEW"

mkdir -p "$OUTDIR"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "old : $OLD"
echo "new : $NEW"
echo "out : $OUTDIR"
echo

## Both sides are flattened. The previous round compared against v9.4, which predates the
## macro system, so only the new side needed expanding. v11.2 carries macros of its own, and
## leaving them unexpanded puts undefined control sequences into the diffed document, which
## is what breaks the citations and mangles the struck-through text.
##
## OLD is expanded against the table set that shipped with it, kept in reviewed/tables_v12_5,
## so a change of table formatting shows up as a change rather than being hidden by both
## sides reading the current file.
OLDTABLES="${OLDTABLES:-$ROOT/reviewed/tables_v12_5}"

flatten() {
  SRC="$1" DST="$2" TABLES="$3" python3 <<'FLAT'
import os, re
src, dst, tables = os.environ['SRC'], os.environ['DST'], os.environ['TABLES']
s = open(src, encoding='utf-8').read()
numbers = None
for tag in re.findall(r'\\input\{(tables/[^}]+)\}', s):
    base = tag.split('/', 1)[1]
    p = os.path.join(tables, base if base.endswith('.tex') else base + '.tex')
    if 'numbers' in tag or 'retention' in tag:
        if os.path.exists(p):
            numbers = p
        s = s.replace('\\input{%s}\n' % tag, '')
        s = s.replace('\\input{%s}' % tag, '')
    else:
        s = s.replace('\\input{%s}' % tag, open(p, encoding='utf-8').read())
macros = {}
for f in sorted(os.listdir(tables)):
    if f.startswith(('numbers', 'retention')) and f.endswith('.tex'):
        for line in open(os.path.join(tables, f), encoding='utf-8'):
            m = re.match(r'\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}\s*$', line)
            if m:
                macros.setdefault(m.group(1), m.group(2))
for name in sorted(macros, key=len, reverse=True):
    s = re.sub(r'\\' + name + r'(?![A-Za-z])', macros[name].replace('\\', '\\\\'), s)
left = sorted({m for m in re.findall(r'\\([A-Za-z]+)', s)} & set(macros))
if left:
    raise SystemExit('unexpanded macros remain in %s: %s' % (os.path.basename(src), left))
if re.search(r'\\input\{', s):
    raise SystemExit('unresolved \\input remains in %s' % os.path.basename(src))
open(dst, 'w', encoding='utf-8').write(s)
print('  flattened %-8s %d macros expanded, no input left' % (os.path.basename(src), len(macros)))
FLAT
}

flatten "$OLD" "$WORK/OLD.tex" "$OLDTABLES"
flatten "$NEW" "$WORK/NEW.tex" "$ROOT/tables"


echo "  running latexdiff"
latexdiff --type=UNDERLINE "$WORK/OLD.tex" "$WORK/NEW.tex" > "$WORK/MARKUP.tex"

echo "  diffing the abstract separately, since the class declares it in the preamble"
WORK="$WORK" python3 <<'PY'
import os, re, subprocess
work = os.environ['WORK']
def abstract_of(p):
    s = open(p, encoding='utf-8').read()
    m = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', s, re.S)
    if not m:
        raise SystemExit('no abstract found in ' + p)
    return m.group(1).strip()

for name, src in [('a_old.tex', 'OLD.tex'), ('a_new.tex', 'NEW.tex')]:
    open(os.path.join(work, name), 'w', encoding='utf-8').write(
        '\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n'
        + abstract_of(os.path.join(work, src)) + '\n\\end{document}\n')

out = subprocess.run(['latexdiff', '--type=UNDERLINE',
                      os.path.join(work, 'a_old.tex'), os.path.join(work, 'a_new.tex')],
                     capture_output=True, text=True)
if out.returncode:
    raise SystemExit(out.stderr)
body = re.search(r'\\begin\{document\}(.*?)\\end\{document\}', out.stdout, re.S).group(1).strip()

## An abstract with no change produces no inline markup, which is the expected outcome in a
## round that did not touch it. Earlier rounds did touch it, and then an empty diff would
## have meant the splice had silently failed, so the two cases are distinguished here rather
## than treated alike.
if 'DIFadd' not in body and 'DIFdel' not in body:
    if abstract_of(os.path.join(work, 'OLD.tex')) != abstract_of(os.path.join(work, 'NEW.tex')):
        raise SystemExit('the abstract differs but the diff produced no inline markup')
    print('  abstract: unchanged this round, nothing to splice')
    raise SystemExit(0)

p = os.path.join(work, 'MARKUP.tex')
mk = open(p, encoding='utf-8').read()
mk2 = re.sub(r'(\\begin\{abstract\}).*?(\\end\{abstract\})',
             lambda m: m.group(1) + '\n' + body + '\n' + m.group(2), mk, count=1, flags=re.S)
if mk2 == mk:
    raise SystemExit('could not splice the abstract markup')
open(p, 'w', encoding='utf-8').write(mk2)
print('  abstract: %d additions, %d deletions marked' % (body.count('DIFadd'), body.count('DIFdel')))
PY

echo "  fitting the wide table and relaxing line breaking"
WORK="$WORK" python3 <<'FITPY'
import os, re
work = os.environ['WORK']
p = os.path.join(work, 'MARKUP.tex')
s = open(p, encoding='utf-8').read()

## A diff carries the old value and the new one in every cell, so the benchmark table is
## about half as wide again as it is in the manuscript and runs past the text block. This
## is a courtesy document rather than a submission, so the table is set smaller and scaled
## to the text width, which keeps every struck and added value legible without reflowing.
m = re.search(r'\\begin\{table\*\}.*?\\end\{table\*\}', s, re.S)
if not m:
    raise SystemExit('no table* found to fit')
blk = m.group(0)
t0 = blk.index(r'\begin{tabular}')
t1 = blk.index(r'\end{tabular}') + len(r'\end{tabular}')
fitted = (blk[:t0]
          + '\\scriptsize\n\\setlength{\\tabcolsep}{3pt}%\n'
          + '\\resizebox{\\textwidth}{!}{%\n' + blk[t0:t1] + '}'
          + blk[t1:])
s = s[:m.start()] + fitted + s[m.end():]

## The same doubling makes several math-heavy lines unbreakable at the manuscript's
## settings, so TeX is given more latitude here than a submission would get.
s = s.replace('\\begin{document}',
              '\\emergencystretch=4em\n\\sloppy\n\\begin{document}', 1)
open(p, 'w', encoding='utf-8').write(s)
print('  table* scaled to the text width, line breaking relaxed')
FITPY

cp "$WORK/MARKUP.tex" "$OUTDIR/$STEM.tex"
for f in wlscirep.cls jabbrv.sty jabbrv-ltwa-all.ldf jabbrv-ltwa-en.ldf naturemag-doi.bst heaped_kde_refs.bib; do
  cp "$ROOT/$f" "$OUTDIR/" 2>/dev/null || true
done
rm -rf "$OUTDIR/figures" && cp -r "$ROOT/figures" "$OUTDIR/figures"

echo "  compiling"
## pdflatex returns nonzero whenever the run had any error, so chaining the passes with &&
## stops after the first one, leaving unresolved citations and a stale bibliography. The
## passes are run unconditionally and the error count is reported at the end instead.
( cd "$OUTDIR"
  rm -f "$STEM.bbl" "$STEM.aux"
  TEXINPUTS=".:" pdflatex -interaction=nonstopmode "$STEM.tex" >/dev/null 2>&1 || true
  BIBINPUTS=".:" BSTINPUTS=".:" bibtex "$STEM" >/dev/null 2>&1 || true
  TEXINPUTS=".:" pdflatex -interaction=nonstopmode "$STEM.tex" >/dev/null 2>&1 || true
  TEXINPUTS=".:" pdflatex -interaction=nonstopmode "$STEM.tex" > "$STEM.build.log" 2>&1 || true )

rm -f "$OUTDIR"/*.aux "$OUTDIR"/*.out "$OUTDIR"/*.blg

echo
ERRS=$(grep -c '^!' "$OUTDIR/$STEM.build.log" || true)
UNDEF=$(grep -c 'Citation.*undefined' "$OUTDIR/$STEM.build.log" || true)
echo "errors:   $ERRS"
echo "undefined citations: $UNDEF"
if [ "$ERRS" != "0" ] || [ "$UNDEF" != "0" ]; then
  echo
  echo "The marked-up file did not compile cleanly. Do not send it out in this state."
  echo "First errors:"
  grep -A2 '^!' "$OUTDIR/$STEM.build.log" | head -12
fi
echo "markup:   $OUTDIR/$STEM.pdf"
echo "source:   $OUTDIR/$STEM.tex"
echo
echo "Blue underline is added text, red strikethrough is deleted text."
