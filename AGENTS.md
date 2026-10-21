# AGENTS.md — 给 AI 助手的协作说明

> 本文件面向在你仓库里工作的 AI 编程助手（以及人类协作者）。**改代码前先读完本文件**。
> 更具体的算法/口径细节见 `题目解读.md`、`scorer/SPEC.md`、`rscore/USAGE.md`。

## 1. 这个项目在做什么

课程评分赛：给定 1000 棵**多叉**基因树（每棵 9000 个类群，类群总并集 10000，内部边标 `S`/`D`），
构造一棵 **10000 叶的二叉物种树 S**，最大化

```
score(S) = Σ_{基因树 g} #{四叶树 q : q 被 S 显示，且在 g 中"可靠"}
```

- **可靠**：`q = ab|cd` 在 g 中，路径 `a→b` 和 `c→d` 全为实线（`S` 或叶子边）。
  等价于"四条末端枝全实线"；锚点之间那条边可以是虚边（见 week5.pdf p45 的图）。
- **显示**：g 中存在一条边，一侧恰含 `{a,b}`、另一侧恰含 `{c,d}`。
- **四叶树是无根的**：分数必须与物种树/基因树的根位置无关（根重挂不变性是硬性回归测试）。
- 打分器算出的分数**就是**比赛分数口径（已与暴力实现逐值对拍）。

## 2. 目录与职责

| 路径 | 内容 | 说明 |
|---|---|---|
| `scorer/` | Python：`ref_brute.py`（暴力，绝对可信）、`fast_formula.py`（快速公式）、`difftest.py` | **唯一权威口径**；改 C++ 前先在这里验证 |
| `scorer/SPEC.md` | 公式推导、两处历史 bug 的修正 | 实现前必读 |
| `rscore/solver/solver.cpp` | C++ 核心：`evalF`（打分）+ NNI 局部搜索 | 主产物 |
| `rscore/rbuild.sh` | 一键：共识初始化 → NNI → 打分 | 对外入口 |
| `rscore/aster2/` | aster2 工作副本 + `src/rscore.hpp`（把计分接进其迭代框架） | 仅用于中等规模基准，**不是主路径** |
| `aster1/` | ASTER v1 上游源码（只读参考，不入库，`.gitignore` 已忽略） | 迭代策略参考：多轮放置 + 三分划 DP 重组 |
| `tools/build_consensus.py` | 簇共识初始化 | 秒级 |
| `submissions/` | 提交树 + `SCOREBOARD.md` 分数台账 | 每次改分必须更新台账 |

## 3. 构建与运行

```bash
# C++ 核心
cd rscore/solver
g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src solver.cpp -o solver

# 打分（10000 类群 × 1000 树 ≈ 2.5 分钟, ~0.3 GB）
./solver -s species.tree -i genes.trees -t 32 --sweeps 0

# 建树（共识 + 3 轮 NNI ≈ 40 分钟）
./rbuild.sh rep01.trees.gz submissions/out.tree 32 3
```

性能基线（32 线程）：打分 2.5 min；NNI 一轮 10 min；初始化打分 90 s；内存约 0.3 GB。

## 4. 编码约定

- C++20；`-O3 -march=native -fopenmp`；不要引入第三方依赖。
- 打分函数 `evalF(gene, bitsL, bitsR, scratch)` 是**性能核心**：O(n) 两趟扫描，必须保持。
- 计数一律 `long long`，总分 `__int128`（单棵树上限 ~1.6e15，总和 ~1.6e18）。
- 多线程：外层 OpenMP 并行，每个线程一个 `Scratch`，**不要**在线程间共享缓冲。
- 物种树必须是二叉；**基因树必须保留多叉**（禁止二值化后再打分，会凭空造四叶树）。

## 5. 已知的坑（踩过，别再踩）

1. **不要用 `newick::Newick::readFile()`** 读基因树：它内部 `arbitrarilyResolve()` 会把多叉二值化，
   并且遇到单叉节点直接抛异常。逐行 `Newick(line)` 构造、自己遍历 `children` 才对。
2. `astral.hpp` 的 `huffman()` 在单叉节点上会崩（`rep01` 每 100 棵树有 11 个单叉节点）。
3. 缺类群：某棵树没有某个类群是正常的（每棵缺 10%），打分时该树只是不贡献涉及它的四叶树。
4. **`cum`（子树内实线成对总数）必须对所有孩子累加，不能按父边实线过滤**（历史 bug）。
5. 计分里"显示"的判据是：
   `(c,d ∉ subtree(u)) 或 (a,b ∉ subtree(v))`，其中 `u=LCA(a,b)`、`v=LCA(c,d)`。
   不要用"祖先关系 + 叶子落在某个孩子里"的旧写法（错的）。
6. `scorer/` 的 Python 实现与 C++ 必须逐值一致；改任何一处都要跑对拍。
7. **计分必须是三分划组合（2026-10-21 修复的 P0）**：不能只数 `Σ_v F(两个孩子簇)`——
   那会漏掉"分裂边在 LCA 之下（3-1 分布）"的四叶树，且分数随根重挂变化。
   正确公式：`W(A,B,C)=F(A,B∪C)+F(B,A∪C)+F(C,A∪B)−2(F(A,B)+F(A,C)+F(B,C))`，`score=ΣW/2`。
   回归测试里必须包含根重挂用例。

## 6. 改完必须做的验证

```bash
./verify.sh        # 见下：构建 + Python 对拍 + 提交树合法性检查
```

`verify.sh` 做三件事：
1. 编译 `solver`；
2. 跑 `scorer/difftest.py`（随机树对拍，与**无根边枚举**真值比对，几百组）；
3. 检查 `submissions/*.tree`：10000 叶、**真二叉**（每个内部节点 2 个孩子）、标签 0–9999 齐全唯一；
4. 40 类群真实数据上 C++ 与 Python 打分一致。

如果改了打分逻辑，**额外**在 40 类群真实数据上对比 C++ 与 Python 的分数（`verify.sh` 里已包含小规模版本）。

## 7. 提交规范

- 分支：`feature/xxx`、`fix/xxx`；
- 提交信息：`feat: ...` / `fix: ...` / `perf: ...` / `docs: ...`；
- 任何影响分数的改动，必须同时更新 `submissions/SCOREBOARD.md`（树、分数、命令、日期）；
- 不要提交：`bin/`、`*.log`、`__pycache__/`、临时数据（`.gitignore` 已覆盖）。

## 8. 当前质量瓶颈（改进方向）

分数高度依赖**初始树**（修正后口径）：500 类群对照里，簇共识 + NNI = 1.889e9，
而 aster2 框架放置 + NNI = 2.535e9（+34%）。全量下我们现有的四棵树彼此只差 0.003%，
说明它们处在同一个盆地，NNI 已经榨不出增量。

下一步计划：受限候选的打分插入（骨架 √(n log n) 个类群打分放置；其余类群在簇伙伴附近
10–20 个候选位置里打分插入），再全量 NNI。详见 `rscore/USAGE.md` 第 9 节。
