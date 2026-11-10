# 数据与产物位置

## 已入库（git）
| 路径 | 内容 | 大小 |
|---|---|---|
| `rep01.trees.gz` | 我们 rep01 的 1000 棵基因树 | 22 MB |
| `启发式算法评分赛全体数据/rep{01..20}.trees.gz` | 全部 rep 的基因树 | 433 MB |
| `submissions/` | 提交树（10000 叶） | ~70 KB/棵 |
| `rscore/results/baseline500_g.trees` | 500 类群 × 100 棵树（基准） | 小 |
| `提交树/评分结果/` | 学生树的评分 JSON | 小 |

## 未入库（需重新生成）
| 路径 | 生成方式 |
|---|---|
| `/tmp/rep01.trees` | `zcat rep01.trees.gz > /tmp/rep01.trees` |
| `/tmp/rep01_spliced.trees` | 去掉单叉节点（ASTRAL 需要），脚本见 `08_LESSONS.md` |
| `/tmp/g1000.trees` | 从 rep16 前 30 棵树取 1000 个类群的子集（脚本见下） |
| `/tmp/astral1000.tree`、`/tmp/incr1000.tree` | 之前的对照树（可重跑） |

## 生成 1000 类群子集

```python
import sys, re, random, gzip
sys.path.insert(0, 'scorer')
from ref_brute import parse_newick
from test_utils import restrict_tree, render
lines = [l.strip() for l in gzip.open('启发式算法评分赛全体数据/rep16.trees.gz', 'rt')][:30]
pat = re.compile(r'(?<=[(,])(\d+)(?=[,)])')
union = set().union(*[set(pat.findall(l)) for l in lines])
rng = random.Random(0); taxa = set(rng.sample(sorted(union), 1000))
out = []
for l in lines:
    ch, lab, nm, par, root = parse_newick(l)
    nc, nl, nn, nr = restrict_tree(ch, lab, nm, root, taxa & set(pat.findall(l)))
    out.append(render(nc, nl, nn, nr) + ";")
open('/tmp/g1000.trees', 'w').write("\n".join(out) + "\n")
```

## 评分结果 JSON 结构
`{"score":..., "n_taxa":..., "n_genes":..., "n_scored_genes":..., "n_skipped_genes":...,
  "n_unused_taxa":..., "n_ignored_taxa":..., "runtime_sec":..., "per_gene":[...]}`
