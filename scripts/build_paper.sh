#!/usr/bin/env bash
# build_paper.sh -- build the paper and the release files.
#
#   dist/ErdosSzekeresHighDim_final2.pdf               the paper
#   dist/ErdosSzekeresHighDim_final2-tex.zip           LaTeX source (main.tex, tables, figures as PNG)
#   dist/ErdosSzekeresHighDim_final2-arxiv.tar.gz      the same source as an arXiv tarball
#   preprintErdosSzekeresHighDim/ErdosSzekeresHighDim.pdf   the paper (the only file there)
#   MasterErdosSzekeresHighDim_masterzipfinal2.zip      every tracked file with dist/ and the
#                                                       preprint folder
#
# Needs pdflatex (REVTeX 4.2, orcidlink), python3, git, zip and tar. The tables in
# paper/tables/ are written by paper/make_tables.py, which rechecks every relation and the
# covering condition; the figures in paper/figures/ by figures/make_figures.py.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
for t in pdflatex python3 git zip tar; do
  command -v "$t" >/dev/null || { echo "missing tool: $t" >&2; exit 2; }
done

NAME=ErdosSzekeresHighDim
PRE=preprint$NAME
MASTER=Master${NAME}_masterzipfinal2.zip

latex() { pdflatex -interaction=nonstopmode -halt-on-error "$@" > /dev/null; }
build() { # build <jobname> <input>, run inside paper/
  rm -f "$1.aux" "$1.out" "$1.log" "$1.toc"
  for run in 1 2 3 4; do
    latex -jobname="$1" "$2"
    [ "$run" -ge 3 ] && ! grep -q "Rerun to get" "$1.log" && break
  done
  if grep -E "Overfull|Underfull|undefined|Rerun to get|stuck" "$1.log"; then
    echo "build_paper.sh: fix the warnings above" >&2
    exit 1
  fi
  rm -f "$1.aux" "$1.log" "$1.out" "$1.toc"
}

echo "== tables"
python3 paper/make_tables.py

echo "== paper"
( cd paper && build main main.tex )

rm -rf dist
mkdir -p dist "$PRE"
cp paper/main.pdf "dist/${NAME}_final2.pdf"
cp paper/main.pdf "$PRE/$NAME.pdf"

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

echo "== tex zip and arXiv tarball"
SRC="$STAGE/${NAME}_final2-tex"
mkdir -p "$SRC/figures" "$SRC/tables"
cp paper/main.tex "$SRC/"
cp paper/figures/*.png "$SRC/figures/"
cp paper/tables/*.tex "$SRC/tables/"
( cd "$SRC" && for i in 1 2 3; do latex -jobname=output main.tex; done \
    && ! grep -E "Overfull|Underfull|undefined" output.log && rm -f output.* outputNotes.bib )
( cd "$STAGE" && zip -qrX "${NAME}_final2-tex.zip" "${NAME}_final2-tex" )
mv "$STAGE/${NAME}_final2-tex.zip" dist/
entries="$(zipinfo -1 "dist/${NAME}_final2-tex.zip" | wc -l)"
[ "$entries" -le 512 ] || { echo "build_paper.sh: the tex zip has $entries entries (at most 512)" >&2; exit 1; }
tar -czf "dist/${NAME}_final2-arxiv.tar.gz" -C "$SRC" main.tex tables figures

echo "== $MASTER"
mkdir -p "$STAGE/master/Master$NAME"
git ls-files -z | grep -zvE "^(dist/|$PRE/|$MASTER\$)" \
  | while IFS= read -r -d '' f; do
      [ -e "$f" ] || continue
      mkdir -p "$STAGE/master/Master$NAME/$(dirname "$f")"
      cp -p "$f" "$STAGE/master/Master$NAME/$f"
    done
mkdir -p "$STAGE/master/Master$NAME/dist" "$STAGE/master/Master$NAME/$PRE"
cp dist/* "$STAGE/master/Master$NAME/dist/"
cp "$PRE/$NAME.pdf" "$STAGE/master/Master$NAME/$PRE/"
rm -f "$MASTER"
( cd "$STAGE/master" && zip -qrX "$ROOT/$MASTER" "Master$NAME" )

echo "built: $(ls dist | tr '\n' ' ')$PRE/$NAME.pdf $MASTER"
