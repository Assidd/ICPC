#!/usr/bin/env bash
set -euo pipefail

IN="${1:-in.txt}"
if [[ ! -f "$IN" ]]; then
  echo "Input file '$IN' not found. Put a test into '$IN' or pass a path as arg."
  exit 1
fi

CXX="${CXX:-g++}"

echo "[1/5] Build solver"
if ! $CXX -O3 -pipe -s -std=c++17 -static main.cpp -o solver 2>/dev/null; then
  $CXX -O3 -pipe -s -std=c++17 main.cpp -o solver
fi

echo "[2/5] Run solver"
./solver < "$IN" > out.txt

echo "[3/5] Build checker"
if ! $CXX -O2 -pipe -s helps/checker.cpp -o checker -static 2>/dev/null; then
  $CXX -O2 -pipe -s helps/checker.cpp -o checker
fi

echo "[4/5] Build scorer"
if ! $CXX -O2 -pipe -s helps/scorer.cpp -o scorer -static 2>/dev/null; then
  $CXX -O2 -pipe -s helps/scorer.cpp -o scorer
fi

echo "[5/5] Validate & score"
./checker "$IN" out.txt
./scorer  "$IN" out.txt
