SP="/home/great_czz/project/演化算法助教工作/算法评分赛/提交树/newtree/rep01_20261116_0950.nwk"
M=128
echo "=== strategy m (task = 1 m) ==="
for t in 1 2 4 8 16 32; do
  ${SECTOR_BIN:-/tmp/sector_sweep} --species "$SP" --genes ${SECTOR_GENES:-/tmp/g.trees} --strategy m --threads $t --m-count $M 2>/dev/null
done
echo "=== strategy block B=8 (task = 8 m's, gene-major) ==="
for t in 1 4 8 16 32; do
  ${SECTOR_BIN:-/tmp/sector_sweep} --species "$SP" --genes ${SECTOR_GENES:-/tmp/g.trees} --strategy block --threads $t --block 8 --m-count $M 2>/dev/null
done
echo "=== strategy nested (all threads on one m, per-thread accumulators) ==="
for t in 1 4 8 16 32; do
  ${SECTOR_BIN:-/tmp/sector_sweep} --species "$SP" --genes ${SECTOR_GENES:-/tmp/g.trees} --strategy nested --threads $t --m-count $M 2>/dev/null
done
