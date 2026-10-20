# 可靠四叶最优树（课程评分赛 · 1 号队伍）

给定一批**多叉、带 S/D 边标号**的基因树，构造一棵**二叉物种树**，使它与所有基因树的
**公共可靠四叶树**数量之和最大。

本仓库包含：题目口径、打分器（参考实现 + 快速实现）、建树工具、实验记录与提交候选树。

## 快速开始

```bash
# 1) 编译核心工具（需要 g++ >= 13）
cd rscore/solver && g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src solver.cpp -o solver && cd ../..

# 2) 给一棵树打分（10000 类群 × 1000 棵树约 2.5 分钟）
./rscore/rbuild.sh --score submissions/team01_nni_v2.tree rep01.trees.gz 32
# → score = 632996196599254

# 3) 从基因树建树（共识初始化 + NNI，约 40 分钟）
./rscore/rbuild.sh rep01.trees.gz submissions/my_tree.tree 32 3

# 4) 跑一遍自检（打分器正确性 + 提交树合法性）
./verify.sh
```

## 当前结果

| 树 | 分数 |
|---|---|
| `submissions/team01_consensus_v0.tree`（簇共识基线） | 533,877,922,026,824 |
| **`submissions/team01_nni_v2.tree`（推荐）** | **632,996,196,599,254** |

详见 [submissions/SCOREBOARD.md](submissions/SCOREBOARD.md)。

## 目录结构

```
├── README.md               ← 本文件
├── AGENTS.md               ← 给 AI 助手的协作说明（队友用 AI 编辑时先读这个）
├── CONTRIBUTING.md         ← 协作流程（分支、提交、验证）
├── verify.sh               ← 一键自检
├── 题目解读.md              ← 题目定义、数据画像、口径确认
├── rep01.trees.gz          ← 1 号队伍的输入数据（1000 棵基因树）
├── 比赛说明.png / week5.pdf ← 题目与讲义
├── scorer/                 ← Python 参考实现（暴力 + 快速公式 + 对拍）
├── rscore/                 ← C++ 工具与实验（主战场）
│   ├── rbuild.sh           ← 一键建树/打分
│   ├── solver/             ← 打分 + NNI 局部搜索（核心）
│   ├── aster2/             ← aster2 工作副本（含 rscore.hpp：把计分接进其迭代框架）
│   ├── results/            ← 实验产物（树、日志）
│   └── USAGE.md            ← 工具使用说明
├── tools/                  ← 早期脚本（簇共识、四叶树语义原型）
└── submissions/            ← 提交候选 + 分数台账
```

## 文档索引

| 想了解 | 看 |
|---|---|
| 题目定义、数据长什么样、口径 | [题目解读.md](题目解读.md) |
| 工具怎么用、参数、性能 | [rscore/USAGE.md](rscore/USAGE.md) |
| 打分公式推导与验证 | [scorer/SPEC.md](scorer/SPEC.md) |
| 当前进展与断点 | [rscore/STATUS.md](rscore/STATUS.md) |
| 用 AI 协作时的约定 | [AGENTS.md](AGENTS.md) |
| 分支/提交/评审流程 | [CONTRIBUTING.md](CONTRIBUTING.md) |

## 协作

```bash
git clone <repo-url> && cd <repo>      # 首次
git checkout -b feature/your-topic     # 开新分支
./verify.sh                            # 改完先自检
git add -A && git commit -m "feat(solver): ..." && git push -u origin feature/your-topic
```
合并前请：① `./verify.sh` 通过；② 若分数有变化，更新 `submissions/SCOREBOARD.md`。
