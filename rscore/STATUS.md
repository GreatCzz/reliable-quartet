# 中断点状态（2026-10-20 17:15）

> 机器已准备关机；所有后台任务均已结束，没有残留进程。

## 1. 当前结果（可直接复现）

**全量 10000 类群**，用「簇共识 v0 → NNI 局部搜索（200 棵树）」得到的树
`results/solver_full_v1.tree`：

| 树 | 全部 1000 棵树的分数 |
|---|---|
| 提交基线 `submissions/team01_consensus_v0.tree` | 533,877,922,026,824 |
| `results/solver_full_v1.tree`（NNI 5 轮） | **632,949,050,483,880**（**+18.6%**） |

树已校验：10000 个叶子、标签 0–9999 齐全且唯一、括号平衡、逗号 9999（二叉）。

**500 类群对照实验**（同一套打分器，Python 独立复核）：

| 树 | 分数 |
|---|---|
| 框架基准（aster2 迭代框架跑出） | 588,288,437 |
| 簇共识 v0（同 500 类群） | 612,537,698 |
| 框架基准 + 3 轮 NNI | **1,273,377,908** |

## 2. 文件地图

```
算法评分赛/
├── scorer/                 # 语义与打分参考实现
│   ├── SPEC.md             # 公式推导 + 两处 bug 的修正
│   ├── ref_brute.py        # 暴力参考（≤12 类群，绝对可信）
│   ├── fast_formula.py     # 快速公式（Python，已验证）
│   ├── difftest.py         # 随机对拍
│   └── bench_real.py       # 测速
├── rscore/
│   ├── README.md           # 方案与进展
│   ├── STATUS.md           # ← 本文件
│   ├── aster2/             # 工作副本（含 rscore.hpp + driver 分支 + makefile 目标）
│   │   ├── bin/rscore      # 框架版（中等规模基准用）
│   │   └── src/rscore.hpp  # 计分接入框架的实现（已验证）
│   ├── solver/solver.cpp   # 全量 NNI 求解器（本项目主产物）
│   └── results/            # 所有跑出来的树、日志、基准数据
└── submissions/
    └── team01_consensus_v0.tree   # 当前基线树（旧的提交候选）
```

## 3. 恢复后怎么继续（命令）

```bash
cd 算法评分赛/rscore/solver
# 重新编译（如果二进制丢失）
g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src solver.cpp -o solver

# 给任意树打分（--sweeps 0）
./solver -s ../results/solver_full_v1.tree -i /tmp/rep01.trees -t 32 --sweeps 0

# 继续优化（从当前最好的树出发）
./solver -s ../results/solver_full_v1.tree -i /tmp/rep01.trees -t 32 --sweeps 5 \
         -o ../results/solver_full_v2.tree
```

注意：`/tmp/rep01.trees`（解压后的基因树）和 `/tmp/rep01_200.trees` 在重启后会丢失，
用 `zcat rep01.trees.gz > /tmp/rep01.trees` 重新生成。

## 4. 建议的下一步（优先级从高到低）

1. **把 `solver_full_v1.tree` 复制进 `submissions/` 作为新的提交候选**（当前它比 v0 高出 18.6%）；
2. 用全部 1000 棵树做 NNI（而不是 200 棵）——每轮约 20 分钟，预计还能再涨；
3. 改进 NNI 的邻域：目前只做单边 NNI，可加"子树交换/SPR"；
4. 处理未解析类群的放置策略（现在它们在初始树里位置随机，虽然影响小）；
5. 与框架基准（500 类群）的差距分析：NNI 从中已经拿到 +116%，说明局部搜索方向正确。

## 5. 已知的坑 / 注意事项

- `rscore/solver/solver` 里的 NNI 当前是"首次改进"策略，接受后立即修改树；
- 求解器输出树的根位置是任意的（四叶树与根无关，不影响分数）；
- `solver` 的 `toNewick` 是递归实现，遇到极深（毛毛虫）树有爆栈风险，目前没遇到；
- 打分器 `scorer/fast_formula.py` 与 C++ 版在 1000+ 组随机用例上一致，可信。
