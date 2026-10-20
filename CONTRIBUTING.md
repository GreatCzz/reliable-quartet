# 协作流程

## 0. 环境要求

- Linux（WSL 可用）、g++ ≥ 13、python3 ≥ 3.10（含 numpy、dendropy，用于参考实现）
- 可选：`zcat`（读 `.gz` 输入）

## 1. 首次克隆

```bash
git clone <repo-url>
cd <repo>/rscore/solver
g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src solver.cpp -o solver
./../../verify.sh        # 确认环境 OK
```

## 2. 分支与提交

| 用途 | 分支名 | 示例 |
|---|---|---|
| 新功能 | `feature/<topic>` | `feature/scored-insertion` |
| 修 bug | `fix/<topic>` | `fix/cum-accumulation` |
| 实验（不一定合并） | `exp/<topic>` | `exp/spr-moves` |

提交信息：

```
feat(solver): 加入受限候选的分数插入
fix(scorer): 修正显示判据
perf(evalF): 用 u32 存 cnt 减少内存带宽
docs: 补充工具用法
```

一次提交只做一件事；**改动算法/参数导致分数变化时，同一次提交里更新 `submissions/SCOREBOARD.md`**。

## 3. 验证（合并前必须通过）

```bash
./verify.sh
```

它做四件事：

1. 编译 `solver`；
2. `scorer/difftest.py`：随机二叉/多叉树 + 随机 S/D，与暴力实现逐值对拍；
3. 检查 `submissions/*.tree` 合法性（10000 叶、二叉、标签齐全唯一）；
4. 40 类群真实数据上 C++ 与 Python 打分一致。

> 改打分核心（`evalF` / `fast_formula.py`）时，**两套实现必须同时改并保持逐值一致**，
> 否则以 `ref_brute.py`（暴力）为准。

## 4. 跑实验的约定

- 大实验（全量 1000 棵树）请把日志和产物放到 `rscore/results/`，日志文件名带日期与参数；
- 结果树如果比台账里的好，复制到 `submissions/` 并更新台账；
- 不要提交超过 100 MB 的中间文件；`.gz` 压缩后再提交。

## 5. AI 协作

- 每次让 AI 改代码前，先让它读 `AGENTS.md`（其中列了口径、构建命令、已知坑、验证清单）；
- 让 AI 在提交前跑 `./verify.sh`；
- 如果 AI 改了打分逻辑，要求它给出"与 Python 参考实现的逐值对拍结果"。

## 6. 评审要点

- 分数是否真的提升（用 `rbuild.sh --score` 复核，别只看程序自报）；
- 是否引入二值化（禁止）、是否破坏多叉/缺类群处理；
- 性能是否可接受（打分 ≤ 3 分钟，单轮 NNI ≤ 15 分钟，内存 ≤ 1 GB）；
- 文档/台账是否同步更新。
