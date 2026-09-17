#!/usr/bin/env bash
# make_markup_v10_7.sh
#
# Builds the marked-up manuscript, submission artifact E-3: the reviewed third submission
# (v9.4) against the current revision, with every change shown inline.
#
# You use zsh and this is a bash script, so invoke it as
#
#     bash make_markup_v10_7.sh
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
#   OLD     the reviewed .tex                (default: ../reviewed/heaped_kde_sr_v9_4.tex)
#   NEW     the current .tex                 (default: the paper beside this script's bundle)
#   OUTDIR  where the markup is written      (default: ../markup)

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
OLD="${OLD:-$ROOT/reviewed/heaped_kde_sr_v9_4.tex}"
NEW="${NEW:-$(ls -t "$ROOT"/heaped_kde_sr_v*.tex 2>/dev/null | grep -v v9_4 | head -1)}"
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

cp "$OLD" "$WORK/OLD.tex"

ROOT="$ROOT" NEW="$NEW" WORK="$WORK" python3 <<'PY'
import os, re
root, new, work = os.environ['ROOT'], os.environ['NEW'], os.environ['WORK']
s = open(new, encoding='utf-8').read()

# inline every \input the paper uses, then expand the emitted number macros
for tag in re.findall(r'\\input\{(tables/[^}]+)\}', s):
    p = os.path.join(root, tag if tag.endswith('.tex') else tag + '.tex')
    if 'numbers' in tag:
        s = s.replace('\\input{%s}\n' % tag, '')
    else:
        s = s.replace('\\input{%s}' % tag, open(p, encoding='utf-8').read())

macros = {}
for line in open(os.path.join(root, 'tables/numbers_v10_7.tex'), encoding='utf-8'):
    m = re.match(r'\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}\s*$', line)
    if m:
        macros[m.group(1)] = m.group(2)
for name in sorted(macros, key=len, reverse=True):
    s = re.sub(r'\\' + name + r'(?![A-Za-z])', macros[name].replace('\\', '\\\\'), s)

left = sorted({m for m in re.findall(r'\\([A-Za-z]+)', s)} & set(macros))
if left:
    raise SystemExit('unexpanded macros remain: %s' % left)
if re.search(r'\\input\{', s):
    raise SystemExit('unresolved \\input remains')
open(os.path.join(work, 'NEW.tex'), 'w', encoding='utf-8').write(s)
print('  flattened: %d macros expanded, no \\input left' % len(macros))
PY

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
if 'DIFadd' not in body or 'DIFdel' not in body:
    raise SystemExit('the abstract diff produced no inline markup')

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
( cd "$OUTDIR" \
  && TEXINPUTS=".:" pdflatex -interaction=nonstopmode "$STEM.tex" >/dev/null \
  && bibtex "$STEM" >/dev/null 2>&1 \
  && TEXINPUTS=".:" pdflatex -interaction=nonstopmode "$STEM.tex" >/dev/null \
  && TEXINPUTS=".:" pdflatex -interaction=nonstopmode "$STEM.tex" > "$STEM.build.log" )

rm -f "$OUTDIR"/*.aux "$OUTDIR"/*.out "$OUTDIR"/*.blg

echo
echo "errors:   $(grep -c '^!' "$OUTDIR/$STEM.build.log" || true)"
echo "markup:   $OUTDIR/$STEM.pdf"
echo "source:   $OUTDIR/$STEM.tex"
echo
echo "Blue underline is added text, red strikethrough is deleted text."
