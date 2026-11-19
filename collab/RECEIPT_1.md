# 回执 #1（缺 1：线性项 ↔ K_v 接口对拍）

## 1. real40 拼装对拍结果
- 通过子树：**38 / 38**
- 位置数：**2510**
- 最大残差：**0**
- 唯一差异：我方金标准 `spr44_formula.cpp` 多算 1 个位置/子树
  （有根树根的两条孩子边在无根树上同一条边），是**我方金标准的冗余**，不是对方的问题。

拼装公式（已验证逐值相等）：
    2Δ(X,e) = 2 * ( Σ_{x∈X}[Σ_{v∈P(e0,e)}(T_x(v,dir_c) − T_x(v,dir_a))] + Σ_{v∈P(e0,e)} K_v(X) )

## 2. sector 表口径
- **1 倍四叶树计数**（不是 2 倍）。定义：
      T_x(m,d) = # 可靠四叶树  x a | b c ，其中 a 落在 m 的第 d 个方向
  抽样 25/25 与 Q 定义一致（`Q_gene_tree_4sets({x}, D_d, D_d', D_d'')`）。
- 文件格式：`int64 hdr[2] = (M, nTaxa)`，随后 `[m][d][taxon]` 的 int64 数组。
  - `M` = 内部节点数（deg≥3），按节点编号升序排列
  - `nTaxa` = **类群名整数值最大值 + 1**（不是类群个数！real40 是 9975）
  - `d ∈ {0,1,2}` = 该节点邻接表顺序 = 孩子(解析顺序)后接父节点；根被抑制时用根的另一个孩子顶替
  - `taxon` = 类群名的整数值
- 位置↔增量：叶子 x 从 e0 侧跨过节点 m 进入方向 c 的增量 = `T_x(m,c) − T_x(m,a)`，a = 朝 e0 的方向
- **重要**：T 表是对**当前完整树 S** 算的（x 仍在 X 内），正好是 A_full 口径，可直接用。

## 3. 核心数据结构（C++）
```cpp
struct AdjTree {
    int n = 0, nTaxa = 0;
    std::vector<int> head, nb, tax, tail;
    std::vector<std::vector<int>> deg;               // deg[u] = 邻接顺序（孩子…然后父）
    int add(int t = -1);
    void edge(int u, int v);
    void build();                                    // 建 CSR: head/nb/tail
    int pos(int u, int v) const;                     // u->v 的 CSR 下标
    int degOf(int u) const;
};
// T 表读取（sector_place.cpp）
//   hdr = (M, nTaxa); 文件布局 [m][d][taxon] -> 内存 [m][taxon][d]
std::vector<long long> T((size_t)M * nTaxa * 3, 0);
```
- kv_engine 的 `node_id` 与我方树 id **完全一致**（同为解析先序），无需映射。
- 我方金标准输出：`X叶子列表 \t pred2`（`SKIPBRUTE=1 DUMP=1`），pred2 = 2×比赛分。

## 4. 机器与时间
- 本机：**32 核 / 15 GB RAM**（跑对拍足够）
- 集群：half QOS **256 核 / 240 GB**，节点 bob/carol/dave
- 距最终提交剩余时间：**待指挥确认（我这边看不到赛程）**

## 5. 给对方的接口请求（复述）
请把线性项做成**可注入入参**：每个子树调用时传 `linear_at_node[v]`（长度 = |S'|，1× 四叶数），
引擎只输出 `node_id \t 2Δ`。这样两边零耦合，我也能用 Python 假数据先验证。
