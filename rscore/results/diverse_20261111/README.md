# 多样化放置树存档（2026-11-11）

数据：rep01（10000 类群 × 1000 棵树）。全部用 `placer2`（扇区公式放置）生成。

| 文件 | 生成方式 | 分数（全量 1000 棵树） |
|---|---|---|
| `div_11.tree` | `placer2 --genes 20 --backbone 300 --order rand --seed 11` | 1,441,448,533,998,459 |
| `div_12.tree` | 同上 seed 12 | 1,440,401,681,716,369 |
| `div_13.tree` | 同上 seed 13 | 1,440,530,992,825,743 |
| `div_14.tree` | 同上 seed 14 | 1,440,653,010,155,353 |
| `div_15.tree` | 同上 seed 15 | 1,441,244,632,988,034 |
| `div_16.tree` | 同上 seed 16 | 1,441,108,688,163,091 |
| `div_recomb.tree` | `recombine`（6 棵 div + rbuild，8 棵树子集） | 1,598,565,084,606,859 |

参考：`submissions/team01_rbuild.tree` = 1,598,565,085,124,791（当前最优）。

## 结论（供后续参考）

- 6 棵随机顺序放置树极差仅 0.07% → **全量下放置顺序不产生多样性**；
- DP 重组选中了 rbuild 的结构 → 放置树没有提供互补信息；
- 若要继续，需要"不同且都足够好"的构造方法（例如 ASTRAL），而不是同一方法的多次随机重启。

## 复现命令

```bash
./rscore/solver/placer2 -i /tmp/rep01.trees -o out.tree --genes 20 --backbone 300 --order rand --seed 11
./rscore/score.sh -t out.tree -g rep01.trees.gz -j 32
```
