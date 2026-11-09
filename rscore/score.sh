#!/usr/bin/env bash
# score.sh — 仅评分模式（助教用）：读入一棵估计物种树 + 一个基因树文件，输出 JSON 分数。
#
# 用法:
#   ./rscore/score.sh -t <species.tree[.gz]> -g <genes.trees[.gz]> [-j THREADS] [-q]
#
# 输出: 一行 JSON（score / n_taxa / n_genes / per_gene ...）
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
BIN="$HERE/solver/score"
[[ -x "$BIN" ]] || { echo "Error: $BIN not found; build it with: cd rscore/solver && g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src score.cpp -o score" >&2; exit 2; }

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
    -h|--help)   sed -n '2,10p' "$0"; exit 0 ;;
    *) echo "Error: unknown argument $1" >&2; exit 2 ;;
  esac
done
[[ -n "$TREE" && -n "$GENES" ]] || { echo "Error: -t <species tree> and -g <gene trees> are required" >&2; exit 2; }

args=(-t "$(plain "$TREE")" -g "$(plain "$GENES")")
[[ -n "$THREADS" ]] && args+=(-j "$THREADS")
[[ -n "$QUIET" ]] && args+=("$QUIET")
exec "$BIN" "${args[@]}"
