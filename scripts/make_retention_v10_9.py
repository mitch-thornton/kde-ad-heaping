#!/usr/bin/env python3
"""make_retention_v10_9.py

Measures how much of the reviewed submission survives into the current revision, and emits
the result as LaTeX macros so the response letter quotes computed numbers rather than typed
ones, on the same footing as every experimental value in the paper.

The editor is being asked to accept a revision that is substantially longer than what was
reviewed. The defensible version of that request rests on where the growth sits, so this
script reports retention section by section rather than as one number.

Method. Both sources are reduced to word sequences with floats, the bibliography, citation
and reference commands, math and markup removed, so the comparison is of prose rather than
of LaTeX. The current source reads its numbers from tables/numbers_v10_*.tex, so its macros
are expanded first; otherwise a macro name would be compared against the number it stands
for and would read as a change. Retention is the total size of the matching blocks found by
difflib.SequenceMatcher, which is the longest common subsequence of words, so a sentence
that was reordered still counts as retained while a rewritten one does not.

Usage, from the bundle root:

    python3 scripts/make_retention_v10_9.py

Environment:
    OLD     the reviewed .tex        (default: reviewed/heaped_kde_sr_v9_4.tex)
    NEW     the current .tex         (default: the heaped_kde_sr_v10_*.tex at the root)
    OUTDIR  where to write           (default: tables)
"""

import difflib
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OLD = os.environ.get("OLD", os.path.join(ROOT, "reviewed", "heaped_kde_sr_v9_4.tex"))
NEW = os.environ.get("NEW", "")
OUTDIR = os.environ.get("OUTDIR", os.path.join(ROOT, "tables"))

if not NEW:
    cands = sorted(glob.glob(os.path.join(ROOT, "heaped_kde_sr_v*.tex")))
    cands = [c for c in cands if "v9_4" not in os.path.basename(c)]
    if not cands:
        sys.exit("no heaped_kde_sr_v*.tex found at the bundle root; set NEW")
    NEW = max(cands, key=os.path.getmtime)

BACK_MATTER = (
    "Funding",
    "Data availability",
    "Code availability",
    "Author contributions statement",
    "Additional information",
)

SECTIONS = [
    ("Abstract", "Abstract"),
    ("Introduction", "Introduction"),
    ("The algebraic-diversity view of binned data", "Algebraic-diversity section"),
    ("Results", "Results"),
    ("Discussion", "Discussion"),
    ("Methods", "Methods"),
]


def expand(path):
    """Read a source, inline its table inputs and expand its emitted number macros."""
    s = open(path, encoding="utf-8").read()
    for tag in re.findall(r"\\input\{(tables/[^}]+)\}", s):
        p = os.path.join(ROOT, tag if tag.endswith(".tex") else tag + ".tex")
        if "numbers" in tag or "retention" in tag:
            s = s.replace("\\input{%s}\n" % tag, "")
            continue
        if os.path.exists(p):
            s = s.replace("\\input{%s}" % tag, open(p, encoding="utf-8").read())
    macros = {}
    for f in glob.glob(os.path.join(ROOT, "tables", "numbers_v10_*.tex")):
        for line in open(f, encoding="utf-8"):
            m = re.match(r"\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}\s*$", line)
            if m:
                macros[m.group(1)] = m.group(2)
    for name in sorted(macros, key=len, reverse=True):
        s = re.sub(r"\\" + name + r"(?![A-Za-z])", macros[name].replace("\\", "\\\\"), s)
    return s


def words(block):
    b = re.sub(r"(?m)^%.*$", "", block)
    b = re.sub(r"\\begin\{table\*?\}.*?\\end\{table\*?\}", "", b, flags=re.S)
    b = re.sub(r"\\begin\{figure\*?\}.*?\\end\{figure\*?\}", "", b, flags=re.S)
    b = re.sub(r"\\bibliography\{.*", "", b, flags=re.S)
    b = re.sub(r"\\(cite|ref|eqref|label|input|includegraphics)\{[^}]*\}", " CITE ", b)
    b = re.sub(r"\$[^$]*\$", " MATH ", b)
    b = re.sub(r"\\[a-zA-Z]+\*?", "", b)
    b = re.sub(r"[{}\\&~]", " ", b)
    return [w for w in b.split() if re.search(r"[A-Za-z]", w)]


def by_section(src):
    out = {}
    s = re.sub(r"(?m)^%.*$", "", src)
    m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", s, re.S)
    out["Abstract"] = words(m.group(1)) if m else []
    parts = re.split(r"\\section\*\{([^}]+)\}", s)
    for i in range(1, len(parts), 2):
        name = parts[i].strip()
        if name in BACK_MATTER:
            continue
        out[name] = out.get(name, []) + words(parts[i + 1])
    return out


def kept(a, b):
    sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
    return sum(blk.size for blk in sm.get_matching_blocks())


def whole(src):
    s = re.sub(r"(?m)^%.*$", "", src)
    i = s.index(r"\begin{abstract}")
    for tag in BACK_MATTER:
        j = s.find("\\section*{%s}" % tag)
        if j > i:
            s = s[:j]
            break
    return words(s[i:])


old_src, new_src = expand(OLD), expand(NEW)
ow, nw = whole(old_src), whole(new_src)
O, N = by_section(old_src), by_section(new_src)

total_kept = kept(ow, nw)
rows, added = [], {}
for key, label in SECTIONS:
    a, b = O.get(key, []), N.get(key, [])
    if not a and not b:
        continue
    k = kept(a, b)
    added[key] = len(b) - k
    rows.append((label, len(a), len(b), k))

growth = sum(added.values())
meth_a, meth_b = O.get("Methods", []), N.get("Methods", [])
meth_kept = kept(meth_a, meth_b)
intro_a, intro_b = O.get("Introduction", []), N.get("Introduction", [])

macros = []


def mac(name, value):
    macros.append("\\newcommand{\\%s}{%s}" % (name, value))


mac("retOldWords", "{:,}".format(len(ow)).replace(",", "{,}"))
mac("retNewWords", "{:,}".format(len(nw)).replace(",", "{,}"))
mac("retKeptWords", "{:,}".format(total_kept).replace(",", "{,}"))
mac("retPctKept", "%.0f" % (100.0 * total_kept / len(ow)))
mac("retDeletedWords", "{:,}".format(len(ow) - total_kept).replace(",", "{,}"))
mac("retMethodsOld", str(len(meth_a)))
mac("retMethodsNew", "{:,}".format(len(meth_b)).replace(",", "{,}"))
mac("retMethodsPctKept", "%.0f" % (100.0 * meth_kept / len(meth_a)))
mac("retMethodsShareOfGrowth", "%.0f" % (100.0 * added.get("Methods", 0) / growth))
mac("retIntroOld", str(len(intro_a)))
mac("retIntroNew", str(len(intro_b)))
mac("retIntroPctKept", "%.0f" % (100.0 * kept(intro_a, intro_b) / len(intro_a)))

os.makedirs(OUTDIR, exist_ok=True)
hdr = [
    "%% retention_v10_9.tex, emitted by scripts/make_retention_v10_9.py. Do not edit by hand.",
    "%%%% %s against %s" % (os.path.basename(OLD), os.path.basename(NEW)),
]
open(os.path.join(OUTDIR, "retention_v10_9.tex"), "w", encoding="utf-8").write(
    "\n".join(hdr + macros) + "\n"
)

print("retention of the reviewed submission in the current revision")
print("%-30s %7s %7s %7s %9s %8s" % ("section", "v9.4", "now", "kept", "kept of", "added"))
for label, a, b, k in rows:
    print("%-30s %7d %7d %7d %8.0f%% %8d" % (label, a, b, k, 100.0 * k / a if a else 0, b - k))
print("%-30s %7d %7d %7d %8.0f%% %8d" % ("TOTAL", len(ow), len(nw), total_kept,
                                         100.0 * total_kept / len(ow), len(nw) - total_kept))
print()
print("-> %s" % os.path.join(OUTDIR, "retention_v10_9.tex"))
for m in macros:
    print("   " + m)
