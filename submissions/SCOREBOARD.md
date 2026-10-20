# 分数台账

> 口径：10000 类群、1000 棵基因树的"公共可靠四叶树"总数。
> 复核命令：`./rscore/rbuild.sh --score <tree> rep01.trees.gz 32`

| 日期 | 树 | 分数 | 怎么来的 | 备注 |
|---|---|---|---|---|
| 2026-10-20 | `team01_consensus_v0.tree` | 533,877,922,026,824 | `tools/build_consensus.py` | 簇共识基线 |
| 2026-10-20 | `team01_nni_v1.tree` | 632,949,050,483,880 | v0 + NNI(200 棵树, 5 轮) | 首个 NNI 版本 |
| 2026-10-20 | **`team01_nni_v2.tree`** | **632,996,196,599,254** | v1 + NNI(1000 棵树, 3 轮) | **当前最好，推荐提交** |
| 2026-10-20 | `team01_rbuild.tree` | 待填 | `./rbuild.sh rep01.trees.gz ... 32 3` | 工具从零跑通验证 |

## 说明

- 分数由 `rscore/solver` 计算；其打分核心与 Python 参考实现（`scorer/fast_formula.py`）逐值对拍一致，
  后者又与暴力实现（`scorer/ref_brute.py`）一致；
- 每次分数变化请追加一行，并在"怎么来的"里写清命令/参数，保证可复现。
