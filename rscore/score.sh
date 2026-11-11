#!/usr/bin/env bash
# score.sh — 仅评分模式（助教用）：读入一棵估计物种树 + 一个基因树文件，输出 JSON 分数。
#
# 用法:
#   ./rscore/score.sh -t <species.tree[.gz]> -g <genes.trees[.gz]> [-j THREADS] [-q]
#
# 输出: 一行 JSON（score / n_taxa / n_genes / per_gene ...）
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
BIN="$HERE/solver/score_fast"
REF="$HERE/solver/score"
USE_REF=0; VERIFY=0

WORK=$(mktemp -d "${TMPDIR:-/tmp}/score.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

plain() {  # $1=path -> prints plain path
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    local out="$WORK/$(basename "${f%.gz}")"
    zcat "$f" > "$out"
    echo "$out"
  else
    echo "$f"
  fi
}

TREE=""; GENES=""; THREADS=""; QUIET=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--tree)   TREE=${2:?}; shift 2 ;;
    -g|--genes)  GENES=${2:?}; shift 2 ;;
    -j|--threads) THREADS=${2:?}; shift 2 ;;
    -q|--quiet)  QUIET="-q"; shift ;;
    --ref)       USE_REF=1; shift ;;
    --verify)    VERIFY=1; shift ;;
    -h|--help)   sed -n '2,10p' "$0"; exit 0 ;;
    *) echo "Error: unknown argument $1" >&2; exit 2 ;;
  esac
done
[[ -n "$TREE" && -n "$GENES" ]] || { echo "Error: -t <species tree> and -g <gene trees> are required" >&2; exit 2; }

[[ -x "$REF" ]] || { echo "Error: $REF not found; build it with: cd rscore/solver && g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src score.cpp -o score" >&2; exit 2; }
[[ -x "$BIN" ]] || { echo "Error: $BIN not found; build it with: cd rscore/solver && g++ -std=c++20 -O3 -march=native -fopenmp score_fast.cpp -o score_fast" >&2; exit 2; }

args=(-t "$(plain "$TREE")" -g "$(plain "$GENES")")
[[ -n "$THREADS" ]] && args+=(-j "$THREADS")
[[ -n "$QUIET" ]] && args+=("$QUIET")

if [[ $VERIFY -eq 1 ]]; then
  fast=$("$BIN" "${args[@]}")
  ref=$("$REF" "${args[@]}")
  fs=$(python3 -c "import json,sys;print(json.loads(sys.argv[1])['score'])" "$fast")
  rs=$(python3 -c "import json,sys;print(json.loads(sys.argv[1])['score'])" "$ref")
  if [[ "$fs" == "$rs" ]]; then
    echo "$fast"
    echo "[verify] fast == ref = $fs" >&2
    exit 0
  else
    echo "[verify] MISMATCH: fast=$fs ref=$rs" >&2
    exit 3
  fi
fi
if [[ $USE_REF -eq 1 ]]; then exec "$REF" "${args[@]}"; else exec "$BIN" "${args[@]}"; fi
