# 受限 SPR 原型（clone oracle）

## 目标

`solver` 当前的 NNI 只改变局部四分关系。本文档记录第一版受限 SPR 原型：它只验证
一次子树重定位的拓扑和局部评分增量，不接入默认 `solver` 搜索路径。

实现位于 `rscore/solver/spr_oracle.cpp`，编译示例：

```bash
g++ -std=c++20 -O2 -fopenmp -I rscore/aster2/src \
  rscore/solver/spr_oracle.cpp -o /tmp/spr_oracle
/tmp/spr_oracle -s species.tree -i genes.trees -t 8 \
  --genes 4 --max-cuts 16 --targets-per-cut 32
```

## P0 拓扑边界

首版只接受以下 SPR：

```text
gp                         gp
 |                          |
 p        --cut/regraft-->  s
 / \
 r   s                         q
                              / \
                             p   ...
                            / \
                           t   r
```

- `r` 是被搬运的子树，`s` 是 `p` 的另一孩子；
- `gp=parent(p)` 必须存在，暂不剪 root child；
- 目标边 `q->t` 必须在 `r` 外部；
- 排除原位置和 no-op；目标边可以是 root->child；
- 复用 `p` 作为新 connector，保持节点 ID 和 cache 行稳定；
- 暂不处理剪 root child 后的 root 规范化。

## 局部性与 oracle

原树与 SPR 后树的 W 只可能在以下节点变化：

```text
{p} ∪ ancestors(gp -> root) ∪ ancestors(q -> root)
```

`r` 子树内部的 A/B 簇和完整 taxa universe 均不变，因此其 W 不变。oracle 对每个 move：

1. clone species tree；
2. 应用 SPR；
3. 重建 clone 的 `bits/sz`；
4. 对 clone 全部内部节点刷新 cache；
5. 只累加理论 affected 集合中的 cache 差异；
6. 检查未 affected 节点差异为零，并检查局部 delta 等于全 cache delta。

增量使用 `ΣW` 口径，官方分数变化为 `delta/2`；因此合法 delta 必须为偶数，但可以为负值，
因为 oracle 会枚举不一定改进的 move。

## 已完成验证

- `smoke_s40.tree` / 40 棵 gene tree：105 个 move，105 通过；
- `builder_fixed500_nni_20261030.tree` / 4 棵 gene tree：409 个 move，409 通过；
- 测试包含多叉 gene tree、S/D 边和缺失类群；
- 未受影响 cache 行逐 gene 保持不变；局部与全量 delta 逐值一致；
- 当前版本还检查 species/gene taxa 集合必须完全相同、root size 必须覆盖全集、叶 bitset/内部 OR bitset、parent/child、内部节点列表和叶 taxon 范围；非法 taxa 输入会以错误码 2 拒绝；
- 另有独立内存 Python 验证：seed `73198`，8 棵随机 7-taxa species tree，每棵配 3 棵随机 2–4 叉 S/D gene tree（1 棵全集、2 棵各缺 1 taxon）；全枚举 80 cuts、726 target edges、2178 gene×target、10206 个 main cache 行，batch delta、clone 全量 W delta 和 `2×ref_brute` 全部一致；该测试未保存源码，仅作为本轮交叉验证摘要，不替代持久化回归测试。## 第三阶段：扩大 cut 覆盖与多轮 SPR（2026-11-11）

搜索 harness 改为每个固定 cut 只计算一次 generalized batch，然后扫描该 cut 的全部合法
目标边；此前的 target 截断只影响 shortlist，不会减少 batch 本身的目标统计成本。

起点固定为 `builder_fixed500_nni_20261030.tree`，500 taxa × 100 genes，所有 SPR 直接 delta
均使用全部 100 棵 gene tree，并在每轮接受一个严格正增益 move 后运行 NNI 到无改进：

| 阶段 | SPR 后分数 | NNI 后分数 |
|---|---:|---:|
| S0 | 2,357,934,671 | 2,357,934,671 |
| 第 1 轮 | 2,360,001,850 | 2,360,100,004 |
| 第 2 轮 | 2,361,433,614 | 2,361,498,973 |
| 第 3 轮 | 2,369,280,736 | **2,371,060,826** |

相对 S0，第三轮最终提升 `13,126,155`（约 `0.557%`）。第三轮 NNI 在第 8 轮停止。
这些是 500×100 实验结果，不代表正式 10000 taxa 结果，也没有修改正式提交树。

第 3 轮使用 32 个 cut（每个 cut-size 桶 4 个）、每个 cut 扫描全部合法 target，screen/full
均为 100 genes，shortlist 为每桶 top-4 加 global top-16。第 1/2 轮同样采用 full exact
batch；每轮均只接受正 delta，再统一 NNI。

后续若扩展实验，应同时保留 exhaustive target delta 才能计算 sampled-vs-exhaustive
recall；不能把随机 target 截断误称为节省 batch 计算。screen genes 只可用于排序候选，
full exact delta 才能接受 move。

固定 cut 子树后，`spr_oracle.cpp` 将 detached 子树视为大小为 `|R|` 的块，在完整 taxa universe
上对 root-connected 主组件一次计算 `dA/dB/dC`、prefix 和每个目标边的 connector 贡献，公式为：

```text
DeltaSigmaW = detachDelta + base + prefix[q]
              + (dSide[q] - dC[q]) + newNode[q->t]
```

其中 `q->t` 是目标有向边，connector 必须是 `(t,R)`；`prefix[q]` 只包含 q 的严格祖先 correction，q 自身方向由显式的
`dSide(q,t)-dC[q]` 加入。数值上 `N_e=W(t,R,U\(t∪R))` 可利用两侧对称性，但实现方向仍必须依据目标边 `t`：不能把 q 的另一孩子当作整个互补侧，也不能使用与 R 重叠的 `t∪R` 作为一侧。真正的方向错误应在 q 非 root 且 q 外侧非空时检查，逐 cache topology 对拍覆盖这类位置。

detach 后 main root 的 `sz/bits` 只是 `N-|R|` 的主组件表示；所有 Wnode 的补集和 size 仍相对于 immutable 的完整 universe `U,N`。root 的 detach baseline `(A_M,B_M,R)` 的 W 不得强制设为零，root 自身 dC 和随后 q 路径 correction 都必须保留。

- 40 taxa / 4 genes / 全部 5108 个合法 move：5108 passed；
- 500 taxa / 4 genes / 16 cuts × 32 targets：409 passed；
- 500 taxa / 100 genes / 16 cuts × 32 targets：409 passed。

`spr_search.cpp` 是独立实验 harness，不改变默认 solver。一次运行固定输入树和 cuts，流程为：

1. 按子树大小桶随机抽取 cuts（`--cuts-per-bin`，由 `--seed` 固定）；
2. 对每个 cut 只计算一次 generalized batch，并扫描该 cut 的合法 target；
3. screen genes 仅用于排序，按每桶 top-`--per-bin` 加 global top-`--global-top` 形成 shortlist；
4. 用全部 gene trees 对 shortlist 做 exact delta，并逐候选 clone-refresh 对拍；
5. 只接受 full exact delta 严格为正的最佳一步 SPR，然后输出树。

默认参数为 `--screen-genes 8 --cuts-per-bin 8 --max-cuts 0 --max-targets 0 --per-bin 4 --global-top 16`。
`--max-cuts 0` 表示不限制分层抽取后的 cut 总数；正数会用相同 seed 从分层 cut 候选中再抽样限额。
`--max-targets 0` 保留每个 cut 的全部合法 target；设为正数只会在 batch 完成后截断候选列表，当前实现不会减少 generalized batch 的统计成本。`--per-bin 0` 可关闭每桶 shortlist，只保留 global top。
screen gene 当前取输入文件前 K 棵，只能筛选候选，不能作为接受依据；正式接受始终使用全部 gene trees。

oracle 与 search 的参数不同，典型用法如下：

```bash
# oracle：验证受影响 cache 行和 generalized batch
spr_oracle -s species.tree -i genes.trees -t 8 \
  --genes 40 --max-cuts 0 --targets-per-cut 0

# search：一次固定-cut 批量扫描并最多应用一步 SPR
spr_search -s species.tree -i genes.trees -t 8 -o spr.tree \
  --screen-genes 100 --cuts-per-bin 4 --max-targets 0 \
  --per-bin 4 --global-top 16 --seed 20261111
```

`spr_search` 本身不运行 NNI 或第二轮。多轮实验应在每轮 SPR 后显式运行 NNI，再以 NNI 输出树重新建立 cache、重新分层选择 cuts 和重新运行 search；不能跨轮复用旧的 `(cutParent,cutChild)` ID。
每轮应记录 SPR 前、SPR 后和 NNI 后分数，并以 full exact score 独立复核；screen 正增益而 full 负增益的候选必须拒绝。
目标随机截断若启用，还应对 exhaustive target 统计 top-k/最佳正增益 recall，不能把当前截断称为 batch 节省。

指定起点 `builder_fixed500_nni_20261030.tree` 的独立 500×100 基线为 **2,357,934,671**。
按固定 seed 依次运行三轮「一次 SPR → NNI 至收敛」，每轮 NNI 后重新抽取 cuts，结果为：

| 阶段 | SPR 后分数 | NNI 后分数 |
|---|---:|---:|
| S0 | 2,357,934,671 | 2,357,934,671 |
| 第 1 轮 | 2,360,001,850 | 2,360,100,004 |
| 第 2 轮 | 2,361,433,614 | 2,361,498,973 |
| 第 3 轮 | 2,369,280,736 | **2,371,060,826** |

三轮使用 seeds `20261111/20261112/20261113`、每桶 4 cuts、全部合法 targets、100 棵
screen/full genes、每桶 top-4 加 global top-16。最终树在同一命令中由独立 `solver --sweeps 0`
重评分，仍为 **2,371,060,826**；`/tmp/pass3_nni.tree` 与独立重评分输出 hash 相同：
`5a31e9ce7d0f4d86bd10d9dcffcc1736ead51ab18ae2f7777c94345e930d0810`。
这只是 500×100 实验，不代表正式 10000 taxa 结果，也没有修改正式提交树。

固定 cut 子树后，`spr_oracle.cpp` 将 detached 子树视为大小为 `|R|` 的块，在完整 taxa
universe 上对 root-connected 主组件计算 `dA/dB/dC` 和 prefix，再扫描目标边的 `newNode`，
公式为：

```text
DeltaSigmaW = detachDelta + base + prefix[t] + newNode[t]
```

`prefix[t]` 已包含到目标父节点 `q` 的严格祖先 correction 以及 q 对 t 方向的
`dSide(q,t)-dC[q]`；connector 必须是 `(t,R)`，不能把 q 的另一孩子当作整个互补侧。
当前 detached main 采用与旧 row 相同的完整-universe W 计算，因此 `dC/base` 数值上恒为
零；代码保留它们是为了使公式与未来增量实现清晰，`dA/dB` 仍用于 prefix。所有 Wnode
的补集和 size 仍相对于 immutable 的完整 universe `U,N`；detached main root 的 `sz/bits`
只是 `N-|R|` 的主组件表示，root 的 detach baseline 和 q 路径 correction 都必须保留。

验证包括 smoke_s40/smoke_g40 全部 5108 个合法 move（`passed=5108, failed=0`），以及
500 taxa 的多组 cut/target 对拍。任何生产化改动都必须保持当前默认 builder/solver 行为，
并与 clone 全量评分逐 move 对拍。

## 10000 taxa 首轮对照（2026-11-11）

输入为 `rep01.trees.gz`（1000 genes），同一只读起点 `submissions/team01_consensus_v0.tree`，24 线程，所有产物限定在 `/tmp/spr-bench-20261111`。这是一轮 NNI sweep 对照，不是 40-sweep 收敛比较。

| 路径 | 结果 | 墙钟 | 峰值 RSS |
|---|---:|---:|---:|
| NNI-only，1 sweep | `1,598,513,034,536,168` → `1,599,795,844,457,596`（+`1,282,809,921,428`，`+0.08025%`；2291 moves） | `38:34.35`，其中初始 scoring 500 s、sweep 1810 s | 278280 KB |
| SPR pilot | 没有接受 move，起点分数不变 | `4797.47 s`（79:57）；user 114663.93 s、sys 38.72 s | 406232 KB |

SPR pilot 使用 seed `20261111`、每个 cut-size 桶 1 cut（共 8 cuts）、全部 target、前 8 棵 gene trees screen、global top-1、1000 棵树 full exact 精评一个候选。159412 个合法 targets 中，screen 预测该候选 `+1,264,859,029`，full exact clone delta 为 `-210,951,707,172`，因此严格拒绝。因无 SPR move 被接受，后续重复 NNI 被取消；不能把此轮记作已完成的 SPR+NNI 质量结果。

该 pilot 的 SPR 搜索耗时约为已测 NNI-only 一 sweep 的 `2.08x`，没有产生分数提升。它证明当前 8-gene screen 的 top-1 在全量数据上可能严重误排，不能据此得出 SPR 本身无效；需要更高质量 screen 或更多 full exact shortlist 候选才能检验 SPR 增益。当前比较也只跑了一次、同为一个 NNI sweep，且 NNI-only 的 user/system CPU 未单独记录，结论限于本次墙钟和这组参数。正式提交树和 SCOREBOARD 均未改动。
