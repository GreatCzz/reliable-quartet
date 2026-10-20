#!/usr/bin/env bash
# rbuild — 从基因树直接构建物种树（簇共识初始化 + NNI 局部搜索），并给出分数。
#
# 用法:
#   ./rbuild.sh <genes.trees[.gz]> [out.tree] [threads] [sweeps]
#
# 例:
#   ./rbuild.sh ../rep01.trees.gz ../submissions/team01.tree 32 3
#
# 只给已有的树打分:
#   ./rbuild.sh --score <species.tree> <genes.trees[.gz]> [threads]
#
# 只精修已有的树（跳过共识初始化）:
#   ./rbuild.sh --refine <species.tree> <genes.trees[.gz]> <out.tree> [threads] [sweeps]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
SOLVER="$HERE/solver/solver"
PY=$(command -v python3)

ensure_plain() {  # $1 = path (maybe .gz) -> prints path to plain file
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    local out="/tmp/$(basename "${f%.gz}")"
    [[ -s "$out" ]] || zcat "$f" > "$out"
    echo "$out"
  else
    echo "$f"
  fi
}

mode="build"
case "${1:-}" in
  --score)  mode="score"; shift ;;
  --refine) mode="refine"; shift ;;
  -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
esac

if [[ "$mode" == "score" ]]; then
  TREE=${1:?usage: rbuild.sh --score <species.tree> <genes[.gz]> [threads]}
  GENES=$(ensure_plain "${2:?usage}")
  THREADS=${3:-32}
  echo -n "score = "
  "$SOLVER" -s "$TREE" -i "$GENES" -t "$THREADS" --sweeps 0 2>&1 | grep "initial score" | sed 's/.*= //' | awk '{print $1}'
  exit 0
fi

if [[ "$mode" == "refine" ]]; then
  INIT=${1:?usage: rbuild.sh --refine <init.tree> <genes[.gz]> <out.tree> [threads] [sweeps]}
  GENES=$(ensure_plain "${2:?}")
  OUT=${3:?}
  THREADS=${4:-32}
  SWEEPS=${5:-3}
else
  GENES=$(ensure_plain "${1:?usage: rbuild.sh <genes[.gz]> [out.tree] [threads] [sweeps]}")
  OUT=${2:-./species.tree}
  THREADS=${3:-32}
  SWEEPS=${4:-3}
  INIT=/tmp/rbuild_init.tree
  echo "[rbuild] 1/3 簇共识初始化 ..."
  "$PY" "$HERE/../tools/build_consensus.py" "$GENES" "$INIT" >/dev/null
fi

echo "[rbuild] 2/3 NNI 局部搜索（$SWEEPS 轮, $THREADS 线程）..."
OMP_NUM_THREADS="$THREADS" "$SOLVER" -s "$INIT" -i "$GENES" -t "$THREADS" --sweeps "$SWEEPS" -o "$OUT" 2>&1 | grep -E "initial score|sweep|final score"

echo "[rbuild] 3/3 输出: $OUT"
echo -n "[rbuild] 最终分数 = "
"$SOLVER" -s "$OUT" -i "$GENES" -t "$THREADS" --sweeps 0 2>&1 | grep "initial score" | sed 's/.*= //' | awk '{print $1}'
