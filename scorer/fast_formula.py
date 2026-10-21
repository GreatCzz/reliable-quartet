#!/usr/bin/env python3
"""快速公式版打分器。

记 S 的内部节点 v 的两个孩子簇为 l, r。每个四叶树在 S 的 LCA 处 2-2 分家，
因此

    score = Σ_v F(l_v, r_v)

对某个簇 X，在基因树 g 上沿"实线可达"做聚合：

    C(u)   = Σ_{孩子 c} [边 u→c 实线] · cnt[c]     # 从 u 沿实线可达、且在 X 中的叶子数
    cnt[u] = C(u)（叶子节点取 1/0）
    anc[u] = Σ_{i<j} c_i·c_j                      # 以 u 为锚点的实线成对数（只看孩子）
    cum[u] = anc[u] + Σ_{孩子} cum[c]             # 子树内实线成对总数
    G[u]   = Σ_{w 是 u 的真祖先, 路径全实线} (C(w) − cnt[ch(w→u)])
    H[u]   = cnt[u]·G[u]                          # 以 u 为锚点的成对里、有一片落在 u 子树内的组合数

四叶树 ab|cd（a,b 在 l，c,d 在 r）在 g 中显示的充要条件是
    (c,d ∉ subtree(u)) 或 (a,b ∉ subtree(v))，其中 u=LCA(a,b), v=LCA(c,d)。

于是

    F(l,r) = cum_l(root)·cum_r(root)          # 全部实线成对组合
             − Σ_u anc_l(u)·anc_r(u)          # u = v（必不显示）
             − Σ_v anc_r(v)·H_l(v)            # u 是 v 的真祖先，且 l 的成对有一片落在 subtree(v)
             − Σ_u anc_l(u)·H_r(u)            # 对称
"""
from ref_brute import parse_newick, leaves_under


def cluster_arrays(ch, lab, root, member, m):
    cnt = [0] * m
    anc = [0] * m
    cum = [0] * m
    for u in range(m):                      # 子节点 id < 父节点 id（本解析器保证）
        if not ch[u]:
            cnt[u] = 1 if member[u] else 0
            continue
        C = 0
        cross = 0
        for c in ch[u]:
            v = cnt[c] if lab[c] != 'D' else 0
            cross += C * v
            C += v
        cnt[u] = C
        anc[u] = cross
        cum[u] = cross + sum(cum[c] for c in ch[u])

    # G[v] = Σ_{w 是 v 的真祖先且路径全实线} (C(w) − cnt[ch(w→v)])
    #      = (C(parent) − cnt[v]) + G[parent]（父边为虚边则整段为 0）
    G = [0] * m
    stack = [root]

    def C_of(u):
        return sum((cnt[c] if lab[c] != 'D' else 0) for c in ch[u])

    Ccache = [0] * m
    while stack:
        u = stack.pop()
        cu = C_of(u)
        Ccache[u] = cu
        for c in ch[u]:
            if lab[c] != 'D':
                G[c] = (cu - cnt[c]) + G[u]
            else:
                G[c] = 0
            stack.append(c)
    H = [cnt[u] * G[u] for u in range(m)]
    return cnt, anc, cum, H


def F_xy(gch, glab, groot, memx, memy, m):
    """F(X,Y)：2 片在 X、2 片在 Y 的可靠四叶树数（X、Y 必须不交）。"""
    _, ax, cx, Hx = cluster_arrays(gch, glab, groot, memx, m)
    _, ay, cy, Hy = cluster_arrays(gch, glab, groot, memy, m)
    return (cx[groot] * cy[groot]
            - sum(ax[u] * ay[u] for u in range(m))
            - sum(ay[u] * Hx[u] for u in range(m))
            - sum(ax[u] * Hy[u] for u in range(m)))


def unrooted_score_detailed(species_newick, gene_newick):
    """无根四叶树总分：Σ_v W(A,B,C)/2（A,B 为 v 的两个孩子簇，C 为其余）。"""
    sch, _, snm, _, sroot = parse_newick(species_newick)
    gch, glab, gnm, _, groot = parse_newick(gene_newick)
    m = len(gch)
    gid = {gnm[u]: u for u in range(m) if not gch[u]}
    leafsets = {u: set(leaves_under(sch, snm, u)) for u in range(len(sch))}
    all_leaves = leafsets[sroot]
    per_node = {}
    total = 0
    for v in range(len(sch)):
        if not sch[v] or len(sch[v]) != 2:
            continue
        A = leafsets[sch[v][0]] & set(gid)
        B = leafsets[sch[v][1]] & set(gid)
        C = (all_leaves & set(gid)) - A - B
        def mem(X):
            a = [False] * m
            for t in X:
                a[gid[t]] = True
            return a
        mA, mB, mC = mem(A), mem(B), mem(C)
        mAUB = [x or y for x, y in zip(mA, mB)]
        mBUC = [x or y for x, y in zip(mB, mC)]
        mAUC = [x or y for x, y in zip(mA, mC)]
        W = (F_xy(gch, glab, groot, mA, mBUC, m) + F_xy(gch, glab, groot, mB, mAUC, m)
             + F_xy(gch, glab, groot, mC, mAUB, m)
             - 2 * (F_xy(gch, glab, groot, mA, mB, m)
                    + F_xy(gch, glab, groot, mA, mC, m)
                    + F_xy(gch, glab, groot, mB, mC, m)))
        per_node[v] = W
        total += W
    return total // 2, per_node


def fast_score_detailed(species_newick, gene_newick):
    sch, _, snm, _, _ = parse_newick(species_newick)
    gch, glab, gnm, _, groot = parse_newick(gene_newick)
    m = len(gch)
    gid = {gnm[u]: u for u in range(m) if not gch[u]}
    total = 0
    per_node = {}
    for v in range(len(sch)):
        if not sch[v] or len(sch[v]) != 2:
            continue
        ml = [False] * m
        for t in leaves_under(sch, snm, sch[v][0]):
            if t in gid:
                ml[gid[t]] = True
        mr = [False] * m
        for t in leaves_under(sch, snm, sch[v][1]):
            if t in gid:
                mr[gid[t]] = True
        _, al, cl, Hl = cluster_arrays(gch, glab, groot, ml, m)
        _, ar, cr, Hr = cluster_arrays(gch, glab, groot, mr, m)
        F = (cl[groot] * cr[groot]
             - sum(al[u] * ar[u] for u in range(m))
             - sum(ar[u] * Hl[u] for u in range(m))
             - sum(al[u] * Hr[u] for u in range(m)))
        per_node[v] = F
        total += F
    return total, per_node


def fast_score(species_newick, gene_newick):
    return unrooted_score_detailed(species_newick, gene_newick)[0]
