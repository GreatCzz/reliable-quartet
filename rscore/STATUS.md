# 断点状态（2026-10-21 16:50）

## 已完成
- 计分改为**无根三分划**（`W(A,B,C)=F(A,B∪C)+F(B,A∪C)+F(C,A∪B)−2(F(A,B)+F(A,C)+F(B,C))`，`score=ΣW/2`）；
  - `scorer/ref_brute.py`：真值改为无根逐边枚举；`fast_formula.py`、`solver.cpp`、`rscore.hpp` 全部同步；
  - 对拍：1000+ 组随机（多叉/缺类群/单叉/D 0–50%/根重挂）与真值一致；`verify.sh` 全绿。
- 框架版 `aster2/bin/rscore` 的 `elementScore` 也改为 W（`score_t` 有符号）；同一无根树不同根形式分数一致。
- 全量（10000×1000）修正后分数：
  - consensus_v0 = 1,598,513,034,536,168
  - nni_v1 = 1,598,549,079,897,629 / nni_v2 = 1,598,519,571,771,453
  - **rbuild = 1,598,565,085,124,791（当前最高）**
- `builder.cpp`（打分插入）+ **节点 W 缓存**：
  - 500 类群：构建 **2.7 s**（budget=60, cand=4, genes=4）→ NNI(40) = **2,117,331,551**
  - 对照：簇共识+NNI = 1,888,878,746；框架树+NNI = 2,534,824,702（**已到 84%**）

## 恢复后的优先事项
1. 跑 `./rscore/verify.sh` 确认环境（编译 + 对拍 + 提交树检查）。
2. 用 builder 做 500 类群参数扫描（budget/cand/genes），目标 ≥ 2.3e9；
3. 全量：`builder → solver --sweeps`，预计构建 20–60 分钟、每轮 NNI ~10 分钟（需先按下面 A 优化或直接用现有速度跑）。
4. 优化项 A（融合三色统计，预计 2.5×）：目前 6 次 `evalF` 各扫一遍全树；
   可合并为 2 趟（3 色 cnt/anc/cum/G/H + 3 组 cross 统计），需严格对拍。

## 命令速查
```bash
cd rscore/solver
./builder -i ../results/baseline500_g.trees -o /tmp/b.tree -t 16 --budget 60 --cand 4 --genes 4
./solver -s /tmp/b.tree -i ../results/baseline500_g.trees -t 8 --sweeps 40 -o /tmp/b_nni.tree
```
全量数据 `/tmp/rep01.trees`（重启后会丢）：`zcat rep01.trees.gz > /tmp/rep01.trees`
