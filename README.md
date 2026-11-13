# 可靠四叶最优树（课程评分赛 · 1 号队伍 / 助教组）

给定一批**多叉、带 S/D 边标号**的基因树，构造一棵**二叉物种树**，使它与所有基因树的
**公共可靠四叶树**数量之和最大。

本仓库包含：题目口径、打分器（参考实现 + 快速实现）、建树/搜索工具、评分台账与提交候选树。

## 新会话从这里开始

接手/续作时先读 **[`handoff/00_START_HERE.md`](handoff/00_START_HERE.md)**——任务目标、当前状态、
关键数字、下一步任务和踩过的坑都在那里。

## 当前结果

**就本次比赛数据（rep01，10000 类群 × 1000 棵树）而言，我们的最好提交：**

| 树 | 分数 | 说明 |
|---|---|---|
| **`submissions/team01_sub_1113b.tree`** | **1,740,983,264,996,131** | **当前最好（11-13 14:31）** |
| `submissions/team01_sub_1113.tree` | 1,738,063,330,826,792 | 11-13 00:53 |
| `submissions/team01_sub_1112_1313.tree` | 1,732,600,430,110,212 | 11-12 13:13 |
| `submissions/team01_sub_1112.tree` | 1,697,691,610,284,838 | 11-12 00:40 |
| `submissions/team01_sub_1111.tree` | 1,674,090,151,256,286 | 11-11 夜 |
| `submissions/team01_nni_v3.tree` | 1,602,087,160,774,962 | 增量 NNI（`nni_incr`）精修 |
| `submissions/team01_rbuild.tree` | 1,598,565,085,124,791 | 共识 + 3 轮 NNI（旧实现） |
| `submissions/team01_consensus_v0.tree` | 1,598,513,034,536,168 | 簇共识基线 |

**横向位置**（各 rep 分数被老师调平，故用"相对增益 = 自家分 ÷ 外来树地板 − 1"比较）：
助教组 **+20.85%**，在已交的 17 个组里**排第 2**；最高 17 组 +21.59%（我们差 0.59%）。

全员评分台账：[`提交树/评分结果/ALL_SCORES.md`](提交树/评分结果/ALL_SCORES.md)；
分组名单与提交状态：[`提交树/评分结果/GROUPS.md`](提交树/评分结果/GROUPS.md)；
单页排名 slide：[`提交树/评分结果/ranking_2026-11-11.pptx`](提交树/评分结果/ranking_2026-11-11.pptx)。

## 快速开始

```bash
# 1) 编译核心工具（g++ >= 13）
cd rscore/solver
g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src solver.cpp   -o solver      # 打分 + NNI
g++ -std=c++20 -O3 -march=native -fopenmp nni_incr.cpp                  -o nni_incr    # 增量 NNI（快 ~13x）
g++ -std=c++20 -O3 -march=native -fopenmp score.cpp                     -o score       # 参考打分（慢）
g++ -std=c++20 -O3 -march=native -fopenmp score_fast.cpp                -o score_fast  # 快速打分（~5 秒/全量）
cd ../..

# 2) 给任意一棵树打分（JSON，含逐基因分数）
./rscore/score.sh -t submissions/team01_sub_1113b.tree -g rep01.trees.gz -j 32
# → 1,740,983,264,996,131

# 3) 快速打分（默认后端，约 5 秒/1000 棵树）
./rscore/solver/score_fast -t submissions/team01_sub_1113b.tree -g rep01.trees.gz -j 32 -q

# 4) 全量 NNI 精修（增量实现，约 3 分钟/轮）
./rscore/solver/nni_incr -s submissions/team01_nni_v3.tree -i rep01.trees.gz -t 26 --sweeps 2 -o out.nwk

# 5) 自检（打分器正确性 + 提交树合法性）
./verify.sh
```

## 工具一览

| 工具 | 作用 | 相对速度 |
|---|---|---|
| `rscore/solver/score_fast.cpp` | **快速打分**（导师 13 状态 DP 移植；与参考实现逐值一致） | ~5 秒/1000 棵树 |
| `rscore/solver/score.cpp` | 参考打分（严格、慢） | ~6 分钟/全量 |
| `rscore/solver/nni_incr.cpp` | **增量 13 状态 DP 驱动的 NNI** | ~3 分钟/轮（旧实现 ~38 分钟） |
| `rscore/solver/solver.cpp` | 打分 + NNI（经典实现，逐值对拍基准） | ~38 分钟/轮 |
| `rscore/solver/spr.cpp` | SPR 局部搜索（与整树打分逐值一致） | 首轮 +2.02% |
| `rscore/solver/builder.cpp` | 打分插入构造 | 分钟级 |
| `rscore/solver/recombine.cpp` | 三分划并集 DP 重组 | 秒级 |
| `rscore/aster2/` | ASTER2 工作副本：`rscore.hpp`（O(n) 计分）、`rscore_incr.hpp`（增量计分，仅放置） | 框架放置 1000 类群 ~10 分钟 |
| `aster1/` | ASTER v1（ASTRAL-IV 对照，只读） | 1000 类群 ~30 分钟 |

## 目录结构

```
├── README.md               ← 本文件
├── AGENTS.md               ← 给 AI 助手的协作说明
├── CONTRIBUTING.md         ← 协作流程
├── verify.sh               ← 一键自检
├── handoff/                ← 新会话交接包（入口：00_START_HERE.md）
├── 题目解读.md              ← 题目定义、数据画像、口径确认
├── rep01.trees.gz          ← 1 号队伍输入数据（1000 棵基因树）
├── scorer/                 ← Python 参考实现（暴力 + 快速公式 + 对拍）
├── rscore/                 ← C++ 工具与实验
│   ├── solver/             ← score / score_fast / solver / nni_incr / spr / builder / recombine
│   ├── aster2/             ← ASTER2 副本 + 计分接口（rscore.hpp / rscore_incr.hpp）
│   ├── NNI_INCR.md         ← 增量 NNI 的设计、验证与收敛轨迹
│   ├── FEASIBILITY_incremental.md ← 增量计分可行性报告
│   └── USAGE.md            ← 工具使用说明
├── submissions/            ← 我方提交候选 + 分数台账
└── 提交树/                  ← 各小组提交 + 统一评分（助教视角）
    ├── repNN/              ← 各组提交树
    └── 评分结果/            ← 逐基因分数、全员台账、排名 slide
```

## 文档索引

| 想了解 | 看 |
|---|---|
| 新会话如何接手 | [handoff/00_START_HERE.md](handoff/00_START_HERE.md) |
| 题目定义、数据、口径 | [题目解读.md](题目解读.md) |
| 打分公式推导与验证 | [scorer/SPEC.md](scorer/SPEC.md) |
| 工具用法、参数、性能 | [rscore/USAGE.md](rscore/USAGE.md) |
| 增量 NNI 的设计与实测 | [rscore/NNI_INCR.md](rscore/NNI_INCR.md) |
| 增量计分可行性（13 状态 DP） | [rscore/FEASIBILITY_incremental.md](rscore/FEASIBILITY_incremental.md) |
| 我方分数历史 | [submissions/TEAM_SCORES.md](submissions/TEAM_SCORES.md) |
| 全员评分（助教视角） | [提交树/评分结果/ALL_SCORES.md](提交树/评分结果/ALL_SCORES.md) |
| 用 AI 协作时的约定 | [AGENTS.md](AGENTS.md) |

## 协作

```bash
git clone <repo-url> && cd <repo>
git checkout -b feature/your-topic
./verify.sh
git add -A && git commit -m "feat(solver): ..." && git push -u origin feature/your-topic
```

合并前请：① `./verify.sh` 通过；② 若分数有变化，更新 `submissions/SCOREBOARD.md` 与 `submissions/TEAM_SCORES.md`。
