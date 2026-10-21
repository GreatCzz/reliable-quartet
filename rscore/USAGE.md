# 使用说明（可靠四叶树：打分 + 建树工具）

## 0. 一句话

给一批**多叉、带 S/D 边标号**的基因树，输出一棵**二叉物种树**，并给出它在题目口径下的分数
（所有基因树的"公共可靠四叶树"总数）。

## 1. 目录与文件

```
算法评分赛/
├── rscore/
│   ├── rbuild.sh              ← 一键脚本（推荐入口）
│   ├── solver/solver          ← C++ 核心（打分 + NNI 局部搜索）
│   ├── solver/solver.cpp
│   ├── aster2/bin/rscore      ← 参考实现：把同一计分接进 aster2 迭代框架（中等规模基准用）
│   └── results/               ← 跑出来的树、日志
├── tools/build_consensus.py   ← 簇共识初始化
├── scorer/                    ← Python 参考实现与验证（ref_brute.py / fast_formula.py / difftest.py）
└── submissions/               ← 提交候选
```

## 2. 编译

```bash
cd rscore/solver
g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src solver.cpp -o solver
```

依赖：g++ ≥ 13（与 aster2 相同）；`-fopenmp` 多线程；无其他库依赖。

## 3. 用法

### 3.1 从基因树建树（最常用）

```bash
./rscore/rbuild.sh rep01.trees.gz submissions/my_tree.tree 32 3
#                         基因树(.gz可)  输出树            线程   NNI轮数
```

流程：簇共识初始化（秒级）→ NNI 局部搜索（每轮约 10 分钟）→ 输出树 + 分数。

### 3.2 只给已有的树打分

```bash
./rscore/rbuild.sh --score submissions/my_tree.tree rep01.trees.gz 32
# 输出：score = 632996196599254
```

### 3.3 在已有树上继续精修

```bash
./rscore/rbuild.sh --refine submissions/my_tree.tree rep01.trees.gz submissions/better.tree 32 3
```

### 3.4 直接用核心程序

```bash
./rscore/solver/solver -s init.tree -i genes.trees -t 32 --sweeps 3 -o out.tree
```

| 参数 | 含义 |
|---|---|
| `-s` | 初始物种树（二叉 Newick） |
| `-i` | 基因树文件（每行一棵，多叉、带 S/D） |
| `-t` | 线程数 |
| `--sweeps N` | NNI 轮数（收敛后自动停） |
| `--rounds R --kick K` | 迭代局部搜索：每轮随机扰动 K 次后再 NNI |
| `-o` | 输出树路径 |

## 4. 输入格式

- **基因树**：Newick，一行一棵；内部边标号 `S`（实线）/`D`（虚线）；**叶子边无需标号，默认实线**；
  允许多叉、允许各棵树的类群集合不同（缺失的类群该树不贡献）。
- **物种树**：二叉 Newick；叶子标签集合必须等于基因树标签的**并集**（漏类群会丢分）。

## 5. 输出

- 标准输出：`score = <整数>`（总分数）
- `-o` 指定的文件：物种树（单行 Newick，二叉）

## 6. 性能（本机 32 线程、10000 类群 × 1000 棵树）

| 操作 | 耗时 | 内存 |
|---|---|---|
| 给一棵树打分 | ≈ 5–6 分钟 | ~0.3 GB |
| NNI 一轮（1000 棵树） | ≈ 10 分钟 | ~0.3 GB |
| 从零建树（共识 + 3 轮 NNI + 打分） | ≈ 40 分钟 | ~0.3 GB |

## 7. 当前成绩（rep01，1000 棵树，10000 类群）

| 树 | 分数 |
|---|---|
| `team01_consensus_v0.tree` | 1,598,513,034,536,168 |
| `team01_rbuild.tree`（当前最高） | **1,598,565,085,124,791** |

> 口径已于 2026-10-21 修正为无根三分划（旧分数作废），详见 `scorer/SPEC.md` 与 `submissions/SCOREBOARD.md`。

## 8. 正确性验证

- C++ 打分的核心函数与 Python 参考实现 `scorer/fast_formula.py` 在 1000+ 组随机用例上逐值一致；
- Python 参考实现又与暴力枚举（`scorer/ref_brute.py`）在 1000+ 组随机用例、以及真实数据抽样上一致；
- 框架版 `rscore` 报告的分数与 Python 独立打分一致（40 类群测试：14773 = 14773）。

## 9. 已知限制与下一步

- **初始树质量决定上限**：修正后口径下 500 类群实验里，簇共识 + NNI = 1.889e9，"打分放置 + NNI" = 2.535e9；
- 下一步计划：把"受限候选的打分插入"加进初始化阶段（骨架 ≈ √(n log n) 个类群打分放置，
  其余类群在簇伙伴附近 10–20 个候选位置里打分插入），预计能把全量分数显著推高；
- 目前 `--rounds/--kick`（迭代局部搜索）在 500 类群上提升有限（<1%）。
