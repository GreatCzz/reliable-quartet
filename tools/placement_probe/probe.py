import itertools, random

# ---------- 树工具 ----------
def parse(s):
    s = s.strip().rstrip(';')
    ch, lab = [], []
    def rec(i):
        u = len(ch); ch.append([]); lab.append(''); i += 1
        while True:
            if s[i] == '(':
                c, i = rec(i); ch[u].append(c)
            else:
                j = i
                while j < len(s) and s[j] not in '(),': j += 1
                leaf = len(ch); ch.append([]); lab.append(s[i:j]); ch[u].append(leaf); i = j
            if s[i] == ',': i += 1; continue
            if s[i] == ')': i += 1; break
        j = i
        while j < len(s) and s[j] not in '(),': j += 1
        lab[u] = s[i:j]; i = j
        return u, i
    root, _ = rec(0)
    return ch, lab, root

def leaves_of(ch, lab, u):
    if not ch[u]: return {lab[u]}
    s = set()
    for c in ch[u]: s |= leaves_of(ch, lab, c)
    return s

def all_splits(ch, lab, root):
    out = []
    def dfs(u):
        if not ch[u]: return {lab[u]}
        s = set()
        for c in ch[u]:
            cs = dfs(c); out.append(frozenset(cs)); s |= cs
        return s
    allL = dfs(root)
    return out, allL

def displayed(splits, allL, pairA, pairB):
    for S in splits:
        T = allL - S
        if set(pairA) <= S and set(pairB) <= T: return True
        if set(pairB) <= S and set(pairA) <= T: return True
    return False

def components(ch, lab, code, root):
    """去掉 D 边后的连通分量（叶子集合）"""
    comps = []
    def dfs(u, inComp):
        if not ch[u]:
            inComp.append(lab[u]); return
        for c in ch[u]:
            if code[c] == 'D':
                newComp = []; dfs(c, newComp); comps.append(set(newComp))
            else:
                dfs(c, inComp)
    cur = []; dfs(root, cur); comps.append(set(cur))
    return comps

def solid_pairs(ch, lab, code, root):
    comps = components(ch, lab, code, root)
    pairs = set()
    for C in comps:
        for a, b in itertools.combinations(sorted(C), 2):
            pairs.add(frozenset((a, b)))
    return pairs

def insert_above(ch, lab, root, p, xname):
    ch2 = [list(c) for c in ch]; lab2 = list(lab)
    par = {c: u for u in range(len(ch2)) for c in ch2[u]}
    par0 = par.get(p, -1)
    w = len(ch2); ch2.append([]); lab2.append('')
    leaf = len(ch2); ch2.append([]); lab2.append(xname)
    ch2[w] = [p, leaf]
    if par0 < 0: root2 = w
    else:
        ch2[par0] = [w if c == p else c for c in ch2[par0]]; root2 = root
    return ch2, lab2, root2

# ---------- 暴力 Δ ----------
def delta_brute(sp, p, xname, genes):
    ch, lab, root = sp
    ch2, lab2, root2 = insert_above(ch, lab, root, p, xname)
    splits2, allL2 = all_splits(ch2, lab2, root2)
    olds = sorted(leaves_of(ch, lab, root))
    tot = 0
    for tri in itertools.combinations(olds, 3):
        a, b, c = tri
        for g in genes:
            chg, labg, rootg = g
            spair = solid_pairs(chg, labg, g[3], rootg)
            # 三个配对
            for (u, v, s, t) in [(a,b,c,xname), (a,c,b,xname), (a,xname,b,c)]:
                pass
    return None

# 重新组织：每个基因预计算 solid pairs / splits
def gene_info(g):
    ch, lab, root, code = g
    return solid_pairs(ch, lab, code, root)

def delta_brute2(sp, p, xname, gene_infos):
    ch, lab, root = sp
    ch2, lab2, root2 = insert_above(ch, lab, root, p, xname)
    splits2, allL2 = all_splits(ch2, lab2, root2)
    olds = sorted(leaves_of(ch, lab, root))
    tot = 0
    for a, b, c in itertools.combinations(olds, 3):
        for pair1, pair2 in [((xname,a),(b,c)), ((xname,b),(a,c)), ((xname,c),(a,b))]:
            if not displayed(splits2, allL2, pair1, pair2): continue
            ok = False
            for spair in gene_infos:
                if frozenset(pair1) in spair and frozenset(pair2) in spair: ok = True; break
            if ok: tot += 1
    return tot

def formula_case12(sp, p, xname, gene_infos, gene_splits, allnames):
    ch, lab, root = sp
    A = leaves_of(ch, lab, p)
    B = set(allnames) - A
    tot = 0
    for spair in gene_infos:
        XA = sum(1 for a in A if frozenset((xname,a)) in spair)
        XB = sum(1 for b in B if frozenset((xname,b)) in spair)
        PA = sum(1 for u,v in itertools.combinations(sorted(A),2) if frozenset((u,v)) in spair)
        PB = sum(1 for u,v in itertools.combinations(sorted(B),2) if frozenset((u,v)) in spair)
        tot += XA*PB + XB*PA
    return tot

if __name__ == '__main__':
    # 小例子
    sp = parse('(((a,b),c),(d,e));')
    genes = []
    random.seed(1)
    def rand_gene(names, pD=0.3):
        items = names[:]; random.shuffle(items)
        def rec(xs):
            if len(xs) == 1: return xs[0]
            k = random.choice([2,2,3]); k = min(k, len(xs))
            parts = [xs[i::k] for i in range(k)]; parts = [p for p in parts if p]
            if len(parts) == 1: return parts[0][0]
            return '(' + ','.join(rec(p) for p in parts) + ')' + ('D' if random.random() < pD else 'S')
        return rec(items) + ';'
    names = ['a','b','c','d','e','x']
    lines = [rand_gene(names) for _ in range(3)]
    gs = []
    for ln in lines:
        ch, lab, root = parse(ln)
        code = {i: lab[i] for i in range(len(lab))}
        gs.append((ch, lab, root, code))
    gene_infos = [gene_info(g) for g in gs]
    oldnames = sorted(leaves_of(*sp[:3]))
    print("物种树叶子:", oldnames)
    for g, ln in zip(gene_infos, lines):
        print("  基因树:", ln)
    # 每个节点尝试插入
    ch, lab, root = sp
    rows = []
    for p in range(len(ch)):
        b = delta_brute2(sp, p, 'x', gene_infos)
        f = formula_case12(sp, p, 'x', gene_infos, None, oldnames)
        rows.append((p, b, f))
    print("  p  brute  case1+2  (差值=case3/4)")
    for p, b, f in rows:
        print(f"  {p:2d} {b:6d} {f:8d}  {b-f:6d}")
