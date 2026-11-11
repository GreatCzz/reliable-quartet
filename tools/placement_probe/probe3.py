import itertools
from probe import parse, leaves_of, all_splits, displayed, insert_above, solid_pairs, delta_brute2

def build(sp):
    ch, lab, root = sp
    par = {c: u for u in range(len(ch)) for c in ch[u]}
    node_of = {lab[u]: u for u in range(len(ch)) if not ch[u]}
    return ch, lab, root, par, node_of

def median(ch, par, root, a, b, c):
    def anc(u):
        s = []
        while True:
            s.append(u)
            if u == root: break
            u = par[u]
        return s
    def lca(u, v):
        su = set(anc(u)); x = v
        while x not in su: x = par[x]
        return x
    def depth(x):
        d = 0
        while x != root: x = par[x]; d += 1
        return d
    return max([lca(a,b), lca(a,c), lca(b,c)], key=depth)

def partner_of(sp, m, p, a, b, c):
    """返回：把 x 插到 p 上方时，四叶树 (x,?) 的搭档是 a/b/c 中的谁"""
    ch, lab, root, par, node_of = sp
    node_m = {'a': a, 'b': b, 'c': c}
    # 组件：以 m 为根的三棵子树
    # 把 m 当作根，求每个节点属于哪个"以 m 的哪个孩子为根"的子树
    def root_at_m():
        # 重建 parent 关系（以 m 为根）
        newpar = {}
        order = [m]
        while order:
            u = order.pop()
            for v in (list(ch[u]) + ([par[u]] if u in par else [])):
                if v not in newpar and v != m:
                    newpar[v] = u; order.append(v)
        return newpar
    newpar = root_at_m()
    def top_child(u):
        # u 在"以 m 为根"的树中的第一层祖先
        x = u
        while newpar.get(x, -1) != m and x != m:
            x = newpar[x]
        return x if x != m else None
    if p == m:
        # 依附在 m 上方：取 parent 方向那棵树
        x = par[m]
        while newpar.get(x, -1) != m and x != m:
            x = newpar[x]
        t = x
    else:
        t = top_child(p)
    # t 的子树里的那个三元组叶子
    for name, node in node_m.items():
        # node 是否在 t 的子树中（以 m 为根）
        x = node
        while x != m and x != t:
            x = newpar.get(x, -1)
            if x == -1: break
        if x == t: return name
    return None

def pred(sp, xname, gene_infos):
    built = build(sp)
    ch, lab, root, par, node_of = built
    old = sorted(node_of)
    contrib = {p: 0 for p in range(len(ch))}
    for A, B, C in itertools.combinations(old, 3):
        a, b, c = node_of[A], node_of[B], node_of[C]
        m = median(ch, par, root, a, b, c)
        for p in range(len(ch)):
            who = partner_of(built, m, p, a, b, c)
            if who is None: continue
            partner = {'a': A, 'b': B, 'c': C}[who]
            others = [x for x in (A, B, C) if x != partner]
            for spair in gene_infos:
                if frozenset((xname, partner)) in spair and frozenset(tuple(others)) in spair:
                    contrib[p] += 1
                    break
    return contrib

if __name__ == '__main__':
    import random
    random.seed(1)
    sp = parse('(((a,b),c),(d,e));')
    names = ['a','b','c','d','e','x']
    def rand_gene(names, pD=0.3):
        items = names[:]; random.shuffle(items)
        def rec(xs):
            if len(xs) == 1: return xs[0]
            k = random.choice([2,2,3]); k = min(k, len(xs))
            parts = [xs[i::k] for i in range(k)]; parts = [p for p in parts if p]
            if len(parts) == 1: return parts[0][0]
            return '(' + ','.join(rec(p) for p in parts) + ')' + ('D' if random.random() < pD else 'S')
        return rec(items) + ';'
    lines = [rand_gene(names) for _ in range(3)]
    gs = []
    for ln in lines:
        ch, lab, root = parse(ln)
        gs.append((ch, lab, root, {i: lab[i] for i in range(len(lab))}))
    gene_infos = [solid_pairs(g[0], g[1], g[3], g[2]) for g in gs]
    pr = pred(sp, 'x', gene_infos)
    print("p  brute  median-pred")
    ok = True
    for p in range(len(sp[0])):
        b = delta_brute2(sp, p, 'x', gene_infos)
        print(f"{p:3d} {b:6d} {pr[p]:10d}")
        if b != pr[p]: ok = False
    print("匹配:", ok)
