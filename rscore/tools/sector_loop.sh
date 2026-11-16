set -e
R="/home/great_czz/project/演化算法助教工作/算法评分赛/rscore"
TREE=/tmp/nni500.tree
BEST=2418965337
for r in $(seq 1 40); do
  /tmp/sector_sweep --species $TREE --genes $R/results/baseline500_g.trees --strategy m --threads 8 --m-count 100000 --dump /tmp/T.txt >/dev/null 2>&1
  /tmp/sector_place --species $TREE --T /tmp/T.txt --apply 100 --out /tmp/s_new.tree >/tmp/place.log 2>&1
  nmv=$(grep -c "^   +" /tmp/place.log || true)
  pred=$(grep -o "predicted gain = .*" /tmp/place.log | head -1)
  s=$($R/solver/solver -s /tmp/s_new.tree -i $R/results/baseline500_g.trees -t 32 --sweeps 0 2>&1 | grep -o "initial score = [0-9]*" | grep -o "[0-9]*")
  echo "round $r: moves=$nmv  $pred  -> score=$s  (+$((s-BEST)))"
  if [ "$s" -le "$BEST" ]; then echo "  no improvement, stop"; break; fi
  BEST=$s
  cp /tmp/s_new.tree /tmp/nni500_best.tree
  TREE=/tmp/nni500_best.tree
done
echo "FINAL: $BEST"
