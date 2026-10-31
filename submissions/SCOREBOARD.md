# 分数台账（最终版）

> **口径（2026-10-21 修正）**：四叶树是**无根**的，分数必须与根位置无关。
> 旧版本只数"rooted LCA 处 2-2 分家"，会低估且随根重挂变化——**旧分数全部作废**（见 issue #1）。
>
> 正确公式：`W(A,B,C)=F(A,B∪C)+F(B,A∪C)+F(C,A∪B)−2(F(A,B)+F(A,C)+F(B,C))`，`score=ΣW/2`。
> 复核命令：`./rscore/rbuild.sh --score <tree> rep01.trees.gz 32`（约 5–6 分钟）

## 全量最终结果（10000 类群 × 1000 棵树）

| 树 | 分数 | 说明 |
|---|---|---|
| `team01_consensus_v0.tree` | 1,598,513,034,536,168 | 簇共识基线 |
| `team01_nni_v1.tree` | 1,598,549,079,897,629 | 旧目标下 NNI |
| `team01_nni_v2.tree` | 1,598,519,571,771,453 | 旧目标下 NNI |
| **`team01_rbuild.tree`** | **1,598,565,085,124,791** | **推荐提交**（工具从零跑出） |
| 按频率插入（毛毛虫） | 1,440,305,865,765,178 | 对照，差 10% |

> **最终复核（2026-10-21 23:33）**：`team01_rbuild.tree` 独立重算 = **1,598,565,085,124,791**（与台账一致），
> 输出时间约 9 分钟（32 线程，机器有其它负载）。

四棵共识派生树彼此只差 **0.003%**——同一盆地；NNI 在修正口径下的增益约 +0.03%。

## 中等规模对照（500 类群 × 100 棵树）

| 树 | 分数 |
|---|---|
| 随机树 + NNI | 1,774,771,175 |
| 簇共识 + NNI | 1,888,878,746 |
| builder（顺序插入）+ NNI | 2,318,553,625 ~ 2,413,257,099（随机种子） |
| 8 棵 builder 重组 | 2,413,257,099（= max，无协同） |
| aster2 框架放置 + NNI | **2,534,824,702**（迄今最高） |

## 修复版 builder + NNI 复测（2026-10-30）

代码版本 `3be2b0b`，输入 `rscore/results/baseline500_g.trees`（500 类群 × 100 棵树），
SHA-256：`2871237a4713b42ad55483c3eb854f1b269c8a244daa8d66516b09b97a0ce7e0`。
以 g++ 15.2.0、`-std=c++20 -O3 -march=native -fopenmp` 从源码编译，两个阶段均为 8 线程；
builder 使用频率排序（不传 `--seed`），只用前 4 棵树决策，NNI 使用全部 100 棵树。

| 阶段 | 产物（位于 `rscore/results/`） | 全部 100 棵树分数 | 本次墙钟耗时 |
|---|---|---:|---:|
| 修复版 builder | `builder_fixed500_20261030.tree` | 1,747,614,212 | 9.32 s |
| NNI（最多 40 轮） | `builder_fixed500_nni_20261030.tree` | **2,357,934,671** | 12.91 s |

命令（`$BUILDER`、`$SOLVER` 指向以上述选项新编译的二进制）：

```bash
"$BUILDER" -i rscore/results/baseline500_g.trees -o rscore/results/builder_fixed500_20261030.tree -t 8 --budget 60 --cand 4 --genes 4
"$SOLVER" -s rscore/results/builder_fixed500_20261030.tree -i rscore/results/baseline500_g.trees -t 8 --sweeps 40 -o rscore/results/builder_fixed500_nni_20261030.tree
"$SOLVER" -s rscore/results/builder_fixed500_nni_20261030.tree -i rscore/results/baseline500_g.trees -t 8 --sweeps 0 -o /tmp/builder-fixed500-rechecked.tree
```

第 35 轮仍接受 1 次移动，第 36 轮 `moves=0` 停止，总计接受 4,826 次移动。
逐轮数据见 `rscore/results/builder_fixed500_nni_20261030.csv`。
初始与精修树均通过 500 叶、标签集合及真二叉检查；最终分数经新进程 C++ 完整重算，
以及 Python `fast_score` 对全部 100 棵树求和，三者一致。
最终树 SHA-256：`a1d1c1c75e478f362da0086de7e9b289ab1f90c2335f73000cbd12a17a154c20`。

NNI 相对初始分数提升 **34.9231%**；比旧版同参数历史记录 **2,360,511,229**
低 **2,576,558（0.1092%）**。构造与精修合计 **22.23 s**（不含编译和独立复核），
峰值 RSS 分别为 5,912 / 6,264 KiB。这是单次运行；旧版 NNI 耗时未记录，不能比较端到端速度。
停止仅表示当前 solver 的候选中无严格改进，不代表全局最优；本次未改动 NNI 实现。
本结果属于 500 类群基准，不替换 10000 类群正式提交树及其分数。

## 结论

- 计分口径已修正并严格对拍（暴力 ↔ 公式 ↔ C++）；
- NNI 的停止轮数依初始树而异（本次修复版 builder 为 36 轮），最终分数仍受构造与搜索邻域影响；
- 重组（三分划 DP）**不会比最好的输入差**，但没有多样化输入时等于取 max；
- 能换盆地的框架式放置在中等规模上跑不动（O(m²·n·T)），这是后续要攻的方向。
