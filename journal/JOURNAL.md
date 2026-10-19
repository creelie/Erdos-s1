# Target journal: The Electronic Journal of Combinatorics

| Item | Details |
|---|---|
| Journal | The Electronic Journal of Combinatorics (EJC) |
| Publisher | the journal itself (diamond open access, no fees) |
| ISSN | 1077-8926 |
| Submission | https://www.combinatorics.org/ojs/index.php/eljc/about/submissions |
| Preprints | the journal accepts papers posted as preprints |
| Class file | `dist/ErdosSzekeresHighDim_journal.pdf` is the journal version in the AMS article class (one column, written by `journal/make_journal.py`); accepted papers are typeset in the journal's class |

## Scope fit

EJC publishes research in all areas of combinatorics, including discrete
geometry, and has published work on the Erdős–Szekeres problem and its
variants. The paper gives explicit constructions that raise the lower bounds
for four Erdős–Szekeres numbers in dimensions three and four, with
certificates checked by hand, by programs and by the Lean 4 kernel.

MSC 2020: primary 52C10; secondary 52A37, 05D10.

## Length and format

- Journal version: 20 pages, one column, six 3D figures and thirteen
  tables (seven of them the certificates); built with no
  overfull or underfull boxes, no undefined references and no stuck floats.
- A short italic abstract with no formulas, four keywords, MSC codes and a
  short table of contents.
- Two authors: Deep Bhattacharjee (corresponding author) and Ushashi
  Bhattacharya, with addresses, e-mail addresses and ORCID iDs.
- One file `main.tex` with the references inline, and the tables in
  `tables/`; no `.bib` or `.bbl` file is needed.
  `dist/ErdosSzekeresHighDim_journal-tex.zip` contains these files and the
  figures as PNG.
- Declarations at the end: "Data and code availability" (the repository and
  its Zenodo archive, and the statement on AI assistance from Claude,
  Anthropic), "Ethical declarations" (funding none, competing interests
  none) and "Author contributions". Appendix A gives the detailed proofs of
  the lemmas and Appendix B the certificates.

## Before submitting

- Reread the journal's current instructions on AI assistance; the paper
  states the use of Claude in the data availability section, and
  `AI_PROVENANCE_ErdosSzekeresHighDim.md` gives the details.
- The covering conditions are checked by computer and by the Lean 4 kernel;
  cite the Zenodo DOI of the release so that referees can run
  `scripts/run_all.sh`.
