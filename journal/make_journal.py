"""make_journal.py -- the journal version of the paper, in the AMS article class.

Reads paper/main.tex (REVTeX 4.2), paper/tables/*.tex and paper/figures/*.png and writes
journal/build/main.tex, journal/build/tables/ and journal/build/figures/: one column,
amsart front matter (addresses, subject classes, keywords), full-width floats turned into
ordinary ones and single-column figures scaled to the text width. The text, the theorems,
the tables and the references are those of the REVTeX version.

Usage: python3 journal/make_journal.py      (then pdflatex main.tex three times in journal/build)
"""
import os
import re
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "paper")
OUT = os.path.join(ROOT, "journal", "build")


def between(s, a, b):
    i = s.index(a) + len(a)
    return s[i:s.index(b, i)]


def booktabs(t):
    """REVTeX's ruledtabular becomes a centred tabular with booktabs rules."""
    t = re.sub(r"\\begin\{ruledtabular\}\n(\\begin\{tabular\}\{[^\n]*\})\n", lambda m: "\\centering\n" + m.group(1) + "\n\\toprule\n", t)
    t = t.replace("\\end{tabular}\n\\end{ruledtabular}", "\\bottomrule\n\\end{tabular}")
    return t.replace("\\\\ \\hline", "\\\\ \\midrule")


def main():
    s = open(os.path.join(SRC, "main.tex"), encoding="utf-8").read()
    title = between(s, r"\title{", "}\n")
    abstract = between(s, "\\itshape\n", "\n")
    keywords = between(s, r"\keywords{", "}\n").replace("; ", ", ")
    date = between(s, r"\date{", "}")
    body = s[s.index(r"\tableofcontents") + len(r"\tableofcontents"):s.index(r"\end{document}")]
    macros = s[s.index(r"\theoremstyle{plain}"):s.index(r"\tolerance")]

    # the side-by-side figure pairs of the two-column version become two figures, one above the other
    body = re.sub(r"\\begin\{figure\*\}\[tp\]\n\\centering\n\\begin\{minipage\}\[t\]\{0.485\\textwidth\}\\centering\n(.*?)"
                  r"\\end\{minipage\}\\hfill\n\\begin\{minipage\}\[t\]\{0.485\\textwidth\}\\centering\n(.*?)\\end\{minipage\}\n\\end\{figure\*\}",
                  lambda m: ("\\begin{figure}[htbp]\n\\centering\n" + m.group(1) + "\\end{figure}\n\\begin{figure}[htbp]\n\\centering\n"
                             + m.group(2) + "\\end{figure}").replace(r"[width=\linewidth]", r"[width=0.5\textwidth]"),
                  body, flags=re.S)
    body = body.replace(r"\begin{figure*}", r"\begin{figure}").replace(r"\end{figure*}", r"\end{figure}")
    body = body.replace(r"\includegraphics[width=\columnwidth]", r"\includegraphics[width=0.55\textwidth]")
    # the nested ceiling of Schoenheim's bound is displayed in one column, where inline it spreads the line
    body = re.sub(r"has at least \$(\\lceil\\frac\{n\}\{s\}.*?\\rceil\\rceil)\$ members",
                  lambda m: "has at least\n\\[" + m.group(1) + "\\]\nmembers", body, count=1)
    # the certificate tables go before the references, which end the paper with the addresses
    tables = body[body.index("\\end{thebibliography}") + len("\\end{thebibliography}"):]
    body = body[:body.index("\\end{thebibliography}") + len("\\end{thebibliography}")] + "\n"
    tables = tables.replace("\\clearpage", "").strip()
    i = body.index("%----------------------------------------------------------------------\n\\begin{thebibliography}")
    body = body[:i] + tables + "\n\\clearpage\n\n" + body[i:]
    body = booktabs(body.replace(r"\begin{figure}[tb]", r"\begin{figure}[htb]"))

    head = r"""\documentclass[11pt]{amsart}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{mathptmx}
\usepackage{amsmath,amssymb,amsthm}
\usepackage[hmargin=0.8in,vmargin=0.75in]{geometry}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{xcolor}
\usepackage[colorlinks=true,linkcolor=blue!50!black,citecolor=blue!50!black,urlcolor=blue!50!black]{hyperref}
\usepackage{orcidlink}
\usepackage[protrusion=true,expansion=false]{microtype}
\usepackage{xurl}
\setcounter{tocdepth}{1}
\renewcommand{\topfraction}{.9}\renewcommand{\bottomfraction}{.6}\renewcommand{\textfraction}{.1}
\renewcommand{\floatpagefraction}{.5}\setcounter{topnumber}{3}\setcounter{totalnumber}{4}

""" + macros + r"""\tolerance=2000
\emergencystretch=3em
\hbadness=2000
\hfuzz=1pt

\begin{document}

\title{""" + title + r"""}

\author{Deep Bhattacharjee\,\orcidlink{0000-0003-0466-750X}}
\address{Formerly, Electro-Gravitational Space Propulsion Laboratory (EGSPL), Bhubaneswar, Odisha 751030, India; formerly, University of Calcutta, Senate House, 87/1 College Street, Kolkata 700073, West Bengal, India}
\email{itsdeep@live.com, d.bhattacharjee@erl-forschung.de}
\thanks{Corresponding author: Deep Bhattacharjee.}

\author{Ushashi Bhattacharya\,\orcidlink{0009-0002-2254-3914}}
\address{Formerly, National Taiwan University, Taipei 106319, Taiwan}
\email{rai.bhattacharya.in@gmail.com}

\date{""" + date + r"""}
\subjclass[2020]{Primary 52C10; Secondary 52A37, 05D10}
\keywords{""" + keywords + r"""}

\begin{abstract}
\itshape """ + abstract + r"""
\end{abstract}

\maketitle
\tableofcontents
"""
    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(os.path.join(OUT, "tables"))
    shutil.copytree(os.path.join(SRC, "figures"), os.path.join(OUT, "figures"),
                    ignore=shutil.ignore_patterns("*.pdf"))
    for f in sorted(os.listdir(os.path.join(SRC, "tables"))):
        t = open(os.path.join(SRC, "tables", f), encoding="utf-8").read()
        t = t.replace(r"\begin{table*}[tp]", r"\begin{table}[tp]").replace(r"\end{table*}", r"\end{table}")
        t = booktabs(t.replace(r"\begin{table}[tb]", r"\begin{table}[htb]"))
        open(os.path.join(OUT, "tables", f), "w", encoding="utf-8").write(t)
    open(os.path.join(OUT, "main.tex"), "w", encoding="utf-8").write(head + body + "\\end{document}\n")
    print("wrote", os.path.relpath(OUT, ROOT))


if __name__ == "__main__":
    main()
