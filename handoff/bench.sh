#!/usr/bin/env bash
# 快速基准：编译 + 500/1000 类群打分 + 可选全量
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT/rscore/solver" || exit 1

echo "== 编译 =="
g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src solver.cpp -o solver &&
g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src score.cpp -o score &&
g++ -std=c++20 -O3 -march=native -fopenmp score_fast.cpp -o score_fast &&
g++ -std=c++20 -O3 -march=native -fopenmp incr_dp.cpp -o incr_dp && echo "  编译完成"

show() { python3 -c "import json,sys; d=json.load(sys.stdin); print(f\"{d['score']:,}  ({d['runtime_sec']:.1f}s, {d['n_genes']} genes)\")"; }

echo "== 500 类群 × 100 棵树 =="
printf "  score_fast : "; ./score_fast -t ../results/baseline500.tree -g ../results/baseline500_g.trees -j 8 -q | show
printf "  solver     : "; ./solver -s ../results/baseline500.tree -i ../results/baseline500_g.trees -t 8 --sweeps 0 2>&1 | grep "initial score" | sed 's/.*= //'

if [ -f /tmp/g1000.trees ]; then
  echo "== 1000 类群 × 30 棵树（/tmp/g1000.trees）=="
  for t in /tmp/incr1000.tree /tmp/astral1000.tree; do
    [ -f "$t" ] && printf "  %-24s " "$(basename "$t")" && ./score_fast -t "$t" -g /tmp/g1000.trees -j 8 -q | show
  done
else
  echo "== 1000 类群数据不存在（见 handoff/07_DATA.md 生成）=="
fi

if [ "${FULL:-0}" = "1" ] && [ -f "$ROOT/submissions/team01_rbuild.tree" ]; then
  echo "== 全量（~6 分钟）=="
  zcat "${1:-$ROOT/rep01.trees.gz}" > /tmp/rep01.trees 2>/dev/null || true
  printf "  score_fast : "; ./score_fast -t "$ROOT/submissions/team01_rbuild.tree" -g /tmp/rep01.trees -j 32 -q | show
fi
