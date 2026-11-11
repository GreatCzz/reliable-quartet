"""用"扇区公式"计算所有插入位置的 Δ，并与暴力对拍。
公式（每个基因 g）：
  对每个内部节点 m（非根）：
    A1,A2 = m 的两个孩子子树叶子集；A3 = 其余
    X(S) = |S ∩ C_x|                       C_x = x 的实线连通分量
    P12 = # 实线对 (A1,A2)；P13、P23 用 deg/pairsInside/cross 推出
    向 subtree(A1) 的所有边加 X1*P23；subtree(A2) 加 X2*P13；其余加 X3*P12
"""
import itertools
from probe import parse, leaves_of, solid_pairs, insert_above, all_splits, displayed, delta_brute2

def build(sp):
    ch, lab, root = sp
    par = {c: u for u in range(len(ch)) for c in ch[u]}
    node_of = {lab[u]: u for u in range(len(ch)) if not ch[u]}
    return ch, lab, root, par, node_of

def comp_of_x(gene, xname):
    ch, lab, root, code = gene
    par = {c: u for u in range(len(ch)) for c in ch[u]}
    # 去掉 D 边后的分量
    comps = []
    def dfs(u, cur):
        if not ch[u]: cur.add(lab[u]); return
        for c in ch[u]:
            if code[c] == 'D':
                nc = set(); dfs(c, nc); comps.append(nc)
            else: dfs(c, cur)
    cur = set(); dfs(root, cur); comps.append(cur)
    for C in comps:
        if xname in C: return C, comps
    return set(), comps

def delta_formula(sp, xname, genes):
    ch, lab, root, par, node_of = build(sp)
    n = len(ch)
    leaves = [u for u in range(n) if not ch[u]]
    name_of = {u: lab[u] for u in leaves}
    sub = {}   # node -> set of leaf names
    def dfs_sub(u):
        if not ch[u]: sub[u] = {lab[u]}; return sub[u]
        s = set()
        for c in ch[u]: s |= dfs_sub(c)
        sub[u] = s; return s
    dfs_sub(root)
    total = {u: 0 for u in range(n)}   # 每个节点（作为"边"）的累计加量
    for g in genes:
        Cx, comps = comp_of_x(g, xname)
        spair = solid_pairs(g[0], g[1], g[3], g[2])
        deg = {}
        for C in comps:
            for a in C: deg[a] = len(C) - 1
        def X(S): return len(S & Cx)
        def P(S, T):
            return sum(1 for u,v in itertools.product(sorted(S), sorted(T)) if frozenset((u,v)) in spair)
        # 每个内部节点的 add（用朴素方式：直接给每个边加，便于验证公式）
        for m in range(n):
            if m == root or not ch[m]: continue
            c1, c2 = ch[m]
            A1, A2 = sub[c1], sub[c2]
            A3 = set(name_of.values()) - A1 - A2
            X1, X2, X3 = X(A1), X(A2), X(A3)
            P12, P13, P23 = P(A1,A2), P(A1,A3), P(A2,A3)
            # 给子树里的每条边加
            def add_subtree(u, v):
                total[u] += v
                for c in ch[u]: add_subtree(c, v)
            add_subtree(c1, X1*P23)
            add_subtree(c2, X2*P13)
            # 其余边（不在 subtree(m) 里，但含 m 自己）
            for u in range(n):
                if u == root: continue
                inside = (u != m) and sub[u] <= sub[m]
                if not inside:
                    total[u] += X3*P12
    return {u: total[u] for u in range(n)}

if __name__ == '__main__':
    import random
    random.seed(11)
    sp = parse('(((a,b),c),(d,e));')
    names = ['a','b','c','d','e','x']
    def rand_gene(names, rng, pD=0.3):
        items = names[:]; rng.shuffle(items)
        def rec(xs):
            if len(xs)==1: return xs[0]
            k = rng.choice([2,2,3]); k=min(k,len(xs))
            parts=[xs[i::k] for i in range(k)]; parts=[p for p in parts if p]
            if len(parts)==1: return parts[0][0]
            return '('+','.join(rec(p) for p in parts)+')'+('D' if rng.random()<pD else 'S')
        return rec(items)+';'
    gs = []
    for _ in range(3):
        ln = rand_gene(names, random)
        ch, lab, root = parse(ln)
        gs.append((ch, lab, root, {i: lab[i] for i in range(len(lab))}))
    gene_infos = [solid_pairs(g[0], g[1], g[3], g[2]) for g in gs]
    f = delta_formula(sp, 'x', gs)
    print("p  brute  formula")
    ok = True
    for p in range(len(sp[0])):
        b = delta_brute2(sp, p, 'x', gene_infos)
        print(f"{p:3d} {b:6d} {f[p]:8d}")
        if b != f[p]: ok = False
    print("匹配:", ok)
