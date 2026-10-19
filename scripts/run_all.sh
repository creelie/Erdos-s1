#!/usr/bin/env bash
# run_all.sh -- runs every check of the four point sets of the paper.
#
# Usage: scripts/run_all.sh             (all steps)
#        scripts/run_all.sh --no-lean   (skip the Lean kernel check)
# Environment:
#   CHECKS="python julia c cpp"   the checks to run (default: all four; "" runs none)
#   SKIP_LEAN=1                   same as --no-lean
# Requirements: lean (Lean 4, through elan; checks/lean/lean-toolchain pins the version)
# for the Lean check; python3, julia, a C compiler and a C++17 compiler for the others.
#
#   python  checks/python/check_all.py: every claim of the paper in exact integer arithmetic
#   c       checks/c/check.c: brute force over all k-subsets, ignoring the certificates
#   cpp     checks/cpp/check.cpp: the certificates (relations, Cramer numerators, covering)
#   julia   checks/julia/check.jl: relations, covering and brute force, in BigInt
#   lean    checks/lean/ErdosSzekeresCertificates.lean: general position, relations and
#           covering, proved by the Lean 4 kernel
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(mktemp -d)"
trap 'rm -rf "$BUILD"' EXIT
CC_BIN="${CC:-cc}"
CXX_BIN="${CXX:-c++}"
CHECKS="${CHECKS-python julia c cpp}"
SKIP_LEAN="${SKIP_LEAN:-0}"
[ "${1:-}" = "--no-lean" ] && SKIP_LEAN=1
HD="$ROOT/highdim"

step() { printf '\n=== %s\n' "$1"; }
need() { command -v "$1" >/dev/null 2>&1 || { echo "missing tool: $1" >&2; exit 2; }; }

for c in $CHECKS; do
  case "$c" in
    python) need python3 ;;
    julia) need julia ;;
    c) need "$CC_BIN" ;;
    cpp) need "$CXX_BIN" ;;
    *) echo "unknown check: $c" >&2; exit 2 ;;
  esac
done

if [ "$SKIP_LEAN" != 1 ]; then
  need lean
  step "Lean 4: kernel check of the four certificates"
  ( cd "$ROOT/checks/lean" && lean ErdosSzekeresCertificates.lean ) 2>&1 | tee "$BUILD/lean.log"
  if grep -qE "error|sorry|warning" "$BUILD/lean.log"; then
    echo "Lean check failed" >&2; exit 1
  fi
  n=0
  for set in X4 X3 Y4 Y3; do
    for thm in general_position relations size; do
      grep -q "'${set}_$thm' does not depend on any axioms" "$BUILD/lean.log" \
        || { echo "unexpected axioms for ${set}_$thm" >&2; exit 1; }
      n=$((n + 1))
    done
    grep -q "'${set}_cover' depends on axioms: \[propext\]" "$BUILD/lean.log" \
      || { echo "unexpected axioms for ${set}_cover" >&2; exit 1; }
    n=$((n + 1))
  done
  echo "$n theorems: twelve depend on no axiom, the four covering theorems on propext only"
  if grep -nwE "sorry|admit|native_decide|bv_decide" "$ROOT/checks/lean/ErdosSzekeresCertificates.lean"; then
    echo "sorry, admit, native_decide or bv_decide found" >&2; exit 1
  fi
  if grep -nE "^[[:space:]]*axiom[[:space:]]" "$ROOT/checks/lean/ErdosSzekeresCertificates.lean"; then
    echo "axiom declaration found" >&2; exit 1
  fi
  echo "no sorry, admit, native_decide, bv_decide or axiom declarations"
  step "Lean 4: the file is generated from the certificates"
  python3 -I "$ROOT/checks/lean/make_lean.py" --stdout | cmp - "$ROOT/checks/lean/ErdosSzekeresCertificates.lean"
  echo "checks/lean/ErdosSzekeresCertificates.lean matches highdim/*.json"
fi

for c in $CHECKS; do
  case "$c" in
    python)
      step "Python: every claim of the paper, exact"
      python3 -I "$ROOT/checks/python/check_all.py" ;;
    julia)
      step "Julia"
      julia "$ROOT/checks/julia/check.jl" ;;
    c)
      step "C: brute force over all 8-subsets"
      "$CC_BIN" -O2 -std=c99 -Wall -o "$BUILD/check_c" "$ROOT/checks/c/check.c"
      "$BUILD/check_c" "$HD/witness_r4_11.json" 4 8
      "$BUILD/check_c" "$HD/witness_r3_15.json" 3 8
      "$BUILD/check_c" "$HD/witness_r4_13_k9.json" 4 9
      "$BUILD/check_c" "$HD/witness_r3_18_k9.json" 3 9
      "$BUILD/check_c" "$HD/scheucher_r3_12.json" 3 7 ;;
    cpp)
      step "C++: the certificates"
      "$CXX_BIN" -O2 -std=c++17 -Wall -o "$BUILD/check_cpp" "$ROOT/checks/cpp/check.cpp"
      "$BUILD/check_cpp" "$HD/witness_r4_11.json" "$HD/cert_r4_11.json" 4 8
      "$BUILD/check_cpp" "$HD/witness_r3_15.json" "$HD/cert_r3_15.json" 3 8
      "$BUILD/check_cpp" "$HD/witness_r4_13_k9.json" "$HD/cert_r4_13_k9.json" 4 9
      "$BUILD/check_cpp" "$HD/witness_r3_18_k9.json" "$HD/cert_r3_18_k9.json" 3 9 ;;
  esac
done

printf '\nAll steps passed.\n'
