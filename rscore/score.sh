#!/usr/bin/env bash
# 用法: ./score.sh <species.tree> [genes.trees]
# 输出该物种树在给定基因树上的官方分数（可靠四叶树总数）。
set -euo pipefail
TREE=${1:?usage: score.sh <species.tree> [genes.trees]}
GENES=${2:-/tmp/rep01.trees}
HERE=$(cd "$(dirname "$0")" && pwd)
"$HERE/solver/solver" -s "$TREE" -i "$GENES" -t "${THREADS:-32}" --sweeps 0 2>&1 | grep "initial score" | sed 's/.*= //'
