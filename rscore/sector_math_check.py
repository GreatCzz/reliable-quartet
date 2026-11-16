#!/usr/bin/env python3
"""
Numerical validation of the *sector* (leaf-placement surface) math.

Levels
------
L0  brute force: score(S_pi) for every attachment point pi, by enumerating all
    quartets and using the SPEC definition of "displayed" + "reliable".
L1  T_x(m,d) by enumeration over triples (a,b,c) with a in direction d of m,
    checking reliability of  x a | b c  directly by definition.
L2  T_x(m,d) by the *component formula* (Theorem 1):

        T_g(m,d) = SUM_{q internal} [ |D_d & C_x(q) & A_x| * CPS(q, i_q(x)) ]

        CPS(q,i) = SUM_{j<l, j,l != i} ( Q_d2(j) Q_d3(l) + Q_d2(l) Q_d3(j) )
        Q_d(j)   = |D_d & C_j(q) & A_q|

L3  rerooting (Theorem 2):  Phi(S_pi) - Phi(S_pi0) == G(pi) - G(pi0), where
        G(R) = SUM_m T_x(m, par_R(m))

L4  x-uniformity (Theorem 3): P = |D_d & C_i & A_x| depends on x only through
    (component of q that contains x, solid cluster of x).

All arithmetic is exact integer arithmetic.
Taxon ids and node ids are kept strictly separate (this bit me once).
"""
import itertools
import random
import sys
from collections import deque


# ------------------------------------------------------------------ tree ----
class Tree:
    __slots__ = ("adj", "tax")

    def __init__(self):
        self.adj = []          # adj[u] = [(v, solid?), ...]
        self.tax = {}          # node -> taxon id  (leaves only)

    def new(self, tax=None):
        self.adj.append([])
        if tax is not None:
            self.tax[len(self.adj) - 1] = tax
        return len(self.adj) - 1

    def edge(self, u, v, solid=True):
        self.adj[u].append((v, solid))
        self.adj[v].append((u, solid))

    def nbrs(self, u):
        return [v for v, _ in self.adj[u]]

    def is_leaf(self, u):
        return u in self.tax

    def t2n(self):
        return {t: u for u, t in self.tax.items()}

    def n2t(self):
        return dict(self.tax)


# ------------------------------------------------------------- builders -----
def rand_binary(n, rng):
    """random binary (unrooted) species tree; taxon i is node i"""
    t = Tree()
    for i in range(n):
        t.new(i)
    comps = list(range(n))
    rng.shuffle(comps)
    while len(comps) > 1:
        a = comps.pop(rng.randrange(len(comps)))
        b = comps.pop(rng.randrange(len(comps)))
        u = t.new()
        t.edge(u, a)
        t.edge(u, b)
        comps.append(u)
    # the last created node has degree 2 -> suppress it (unrooted binary tree)
    root = comps[0]
    if not t.is_leaf(root) and len(t.adj[root]) == 2:
        c1, c2 = t.nbrs(root)
        t.adj[c1] = [(v, s) for v, s in t.adj[c1] if v != root]
        t.adj[c2] = [(v, s) for v, s in t.adj[c2] if v != root]
        t.adj[root] = []
        t.edge(c1, c2)
    return t


def rand_gene(taxa, rng, p_dashed=0.30):
    """random *polytomous* gene tree; internal edges solid w.p. 1-p_dashed"""
    children = {}
    nid = 0
    roots = []
    for a in taxa:
        children[nid] = ("leaf", a)
        roots.append(nid)
        nid += 1
    while len(roots) > 1:
        k = min(len(roots), rng.choice([2, 2, 2, 3, 3, 4]))
        grp = [roots.pop(rng.randrange(len(roots))) for _ in range(k)]
        children[nid] = ("int", grp)
        roots.append(nid)
        nid += 1
    root = roots[0]

    t = Tree()
    for _ in range(nid):
        t.new()
    stack = [root]
    while stack:
        u = stack.pop()
        if children[u][0] == "leaf":
            t.tax[u] = children[u][1]
            continue
        for c in children[u][1]:
            solid = True if children[c][0] == "leaf" else (rng.random() > p_dashed)
            t.edge(u, c, solid)
            stack.append(c)
    return t


# --------------------------------------------------------------- helpers ----
def comp_nodes(t, q, w):
    """node set of the component of t-q that contains w"""
    seen = {q}
    st = [w]
    out = set()
    while st:
        u = st.pop()
        if u in seen:
            continue
        seen.add(u)
        out.add(u)
        for v in t.nbrs(u):
            if v not in seen:
                st.append(v)
    return out


def comp_taxa(t, q, w):
    return frozenset(t.tax[u] for u in comp_nodes(t, q, w) if u in t.tax)


def all_dirs(t, q):
    """[(neighbour, taxon set), ...] for the components of t - q"""
    return [(w, comp_taxa(t, q, w)) for w in t.nbrs(q)]


def path_edges(t, a, b):
    """edges (u,v,solid) on the path between NODES a and b"""
    prev = {a: None}
    dq = deque([a])
    while dq:
        u = dq.popleft()
        if u == b:
            break
        for v, s in t.adj[u]:
            if v not in prev:
                prev[v] = (u, s)
                dq.append(v)
    assert b in prev, "disconnected"
    res = []
    cur = b
    while prev[cur] is not None:
        u, s = prev[cur]
        res.append((u, cur, s))
        cur = u
    return res


def topo_S(t, a, b, c, d):
    """NODES a,b,c,d -> the pairing whose two pair-paths are edge disjoint"""
    for (p, q), (r, s) in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))):
        e1 = {frozenset((u, v)) for u, v, _ in path_edges(t, p, q)}
        e2 = {frozenset((u, v)) for u, v, _ in path_edges(t, r, s)}
        if not (e1 & e2):
            return ((p, q), (r, s))
    raise RuntimeError("species tree is not binary / quartet degenerate")


def solid_path(t, a, b):
    return all(s for _, _, s in path_edges(t, a, b))


def displays(t, an, bn, cn, dn):
    """gene tree t displays an bn | cn dn   (NODES)"""
    for u in range(len(t.adj)):
        for v, _ in t.adj[u]:
            if u < v:
                tx = comp_taxa(t, u, v)
                rest = frozenset(t.tax.values()) - tx
                if ({t.tax.get(an), t.tax.get(bn)} <= tx
                        and {t.tax.get(cn), t.tax.get(dn)} <= rest):
                    return True
                if ({t.tax.get(cn), t.tax.get(dn)} <= tx
                        and {t.tax.get(an), t.tax.get(bn)} <= rest):
                    return True
    return False


def reliable(t, an, bn, cn, dn):
    """quartet an bn | cn dn reliable in gene tree t ?  (NODES)"""
    if not displays(t, an, bn, cn, dn):
        return False
    return solid_path(t, an, bn) and solid_path(t, cn, dn)


def solid_cluster(t, x):
    """taxa solid-reachable from NODE x"""
    seen = {x}
    st = [x]
    tx = set()
    while st:
        u = st.pop()
        for v, s in t.adj[u]:
            if s and v not in seen:
                seen.add(v)
                st.append(v)
                if v in t.tax:
                    tx.add(t.tax[v])
    if x in t.tax:
        tx.add(t.tax[x])
    return frozenset(tx)


def solid_visible_from(t, q, w):
    """taxa u in the component of t-q through w with path u->q fully solid"""
    if not any(v == w and s for v, s in t.adj[q]):
        return frozenset()
    seen = {q}
    st = [w]
    tx = set()
    while st:
        u = st.pop()
        if u in seen:
            continue
        seen.add(u)
        if u in t.tax:
            tx.add(t.tax[u])
        for v, s in t.adj[u]:
            if v not in seen and s:
                st.append(v)
    return frozenset(tx)


# ------------------------------------------------------------------- L0 -----
def score_brute(S, genes):
    Sn = S.t2n()
    tot = 0
    for g in genes:
        gn = g.t2n()
        alltax = frozenset(g.tax.values())
        for quad in itertools.combinations(sorted(alltax), 4):
            ((p, q), (r, s)) = topo_S(S, *[Sn[t] for t in quad])
            pairs = [((p, q), (r, s))]
            # map species-node pairing back to taxa
            sp = {u: t for t, u in Sn.items()}
            (p1, q1), (r1, s1) = (sp[p], sp[q]), (sp[r], sp[s])
            if reliable(g, gn[p1], gn[q1], gn[r1], gn[s1]):
                tot += 1
    return tot


# ------------------------------------------------------ species surgery -----
def detach(S, x):
    """S -> S- : remove leaf NODE x, suppress its neighbour"""
    p = S.nbrs(x)[0]
    nb = [w for w in S.nbrs(p) if w != x]
    assert len(nb) == 2, "species tree must be binary"
    nid = {}
    t = Tree()
    for u in range(len(S.adj)):
        if u in (x, p):
            continue
        nid[u] = t.new(S.tax.get(u))
    for u in range(len(S.adj)):
        if u in (x, p):
            continue
        for v, _ in S.adj[u]:
            if u < v and v not in (x, p):
                t.edge(nid[u], nid[v])
    t.edge(nid[nb[0]], nid[nb[1]])
    return t, nid


def attach(Sm, xtax, u, v):
    """Sm with edge (u,v) subdivided; a new leaf with taxon xtax hung on it"""
    t = Tree()
    for i in range(len(Sm.adj)):
        t.new(Sm.tax.get(i))
    lx = t.new(xtax)
    z = t.new()
    for i in range(len(Sm.adj)):
        for w, _ in Sm.adj[i]:
            if i < w and not ((i == u and w == v) or (i == v and w == u)):
                t.edge(i, w)
    t.edge(z, u)
    t.edge(z, v)
    t.edge(z, lx)
    return t


# ---------------------------------------------------------------- L1/L2 -----
def T_direct(g, xnode, Ds, d):
    """#{ (a,b,c) : a in Ds[d], b in Ds[d2], c in Ds[d3],  xa|bc reliable }"""
    d2, d3 = [i for i in range(3) if i != d]
    gn = g.t2n()
    xtax = g.tax[xnode]
    cnt = 0
    for a in Ds[d] & frozenset(gn):
        for b in Ds[d2] & frozenset(gn):
            for c in Ds[d3] & frozenset(gn):
                if reliable(g, xnode, gn[a], gn[b], gn[c]):
                    cnt += 1
    return cnt


def T_formula(g, xnode, Ds, d):
    """Theorem 1 component formula"""
    d2, d3 = [i for i in range(3) if i != d]
    Ax = solid_cluster(g, xnode)
    internal = [q for q in range(len(g.adj)) if not g.is_leaf(q)]
    total = 0
    for q in internal:
        dirs = all_dirs(g, q)
        i0 = None
        for idx, (w, tx) in enumerate(dirs):
            if xnode in g.tax and g.tax[xnode] in tx:
                i0 = idx
                break
        if i0 is None:
            continue
        P = len(Ds[d] & dirs[i0][1] & Ax)
        if P == 0:
            continue
        q2, q3 = [], []
        for idx, (w, tx) in enumerate(dirs):
            if idx == i0:
                continue
            vis = solid_visible_from(g, q, w)
            q2.append(len(Ds[d2] & tx & vis))
            q3.append(len(Ds[d3] & tx & vis))
        cps = 0
        for j in range(len(q2)):
            for l in range(j + 1, len(q2)):
                cps += q2[j] * q3[l] + q3[j] * q2[l]
        total += P * cps
    return total


# ------------------------------------------------------------- L4 check -----
def check_uniformity(g, xnode, Ds):
    """P must depend on x only through (component, solid cluster)"""
    for q in (q for q in range(len(g.adj)) if not g.is_leaf(q)):
        dirs = all_dirs(g, q)
        for idx, (w, tx) in enumerate(dirs):
            bucket = {}
            for y in tx:
                yn = g.t2n()[y]
                cl = solid_cluster(g, yn)
                bucket.setdefault(cl, set()).add(len(tx & cl))
            for cl, vs in bucket.items():
                if len(vs) > 1:
                    return f"P not uniform: q={q} comp={idx} cluster={sorted(cl)}"
    return None


# ------------------------------------------------------------------ L5 ------
def sweep_T(g, Ds, d, gamma, root):
    """Theorem 3: for every x in `gamma`,  T_x(m,d) = sum over the root->x path of
    TWO(u) = EDG(par->u) + EDG(u->par).   Returns {taxon: T}, cluster-restricted."""
    d2, d3 = [i for i in range(3) if i != d]
    internal = [q for q in range(len(g.adj)) if not g.is_leaf(q)]
    dirs_of = {q: all_dirs(g, q) for q in internal}
    EDG = {}
    for q in internal:
        dirs = dirs_of[q]
        q2, q3 = [], []
        for w, tx in dirs:
            vis = solid_visible_from(g, q, w)
            q2.append(len(Ds[d2] & tx & vis))
            q3.append(len(Ds[d3] & tx & vis))
        for i0, (w, tx) in enumerate(dirs):
            P = len(Ds[d] & tx & gamma)
            if P == 0:
                EDG[(q, i0)] = 0
                continue
            cps = 0
            for j in range(len(dirs)):
                if j == i0:
                    continue
                for l in range(j + 1, len(dirs)):
                    if l == i0:
                        continue
                    cps += q2[j] * q3[l] + q3[j] * q2[l]
            EDG[(q, i0)] = P * cps

    def edg(q, w):
        if g.is_leaf(q):
            return 0
        for idx, (ww, _) in enumerate(dirs_of[q]):
            if ww == w:
                return EDG[(q, idx)]
        raise KeyError

    par = {root: None}
    order = [root]
    st = [root]
    while st:
        u = st.pop()
        for v in g.nbrs(u):
            if v not in par:
                par[v] = u
                order.append(v)
                st.append(v)
    # T_x = SUM_q EDG(q, dir_q(x))  splits into
    #   (a) q on the root->x path : dir = child towards x
    #   (b) q not on it (x not in subtree(q)) : dir = parent
    # so  T_x = TOTALP + SUM over path edges (p->u) of [ EDG(p->u) - GP(p) ]
    GP = {q: edg(q, par[q]) for q in order if q != root and not g.is_leaf(q)}
    TOTALP = sum(GP.values())
    acc = {root: TOTALP}
    for u in order:
        if u == root:
            continue
        p = par[u]
        acc[u] = acc[p] + edg(p, u) - GP.get(p, 0)
    return {g.tax[x]: acc[x] for x in g.tax if g.tax[x] in gamma}


# ----------------------------------------------------------- check driver ---
def check_instance(seed, n=8, ngenes=3, p_dashed=0.30, p_missing=0.10):
    rng = random.Random(seed)
    S = rand_binary(n, rng)
    genes = []
    for _ in range(ngenes):
        taxa = [a for a in range(n) if rng.random() > p_missing]
        if len(taxa) < 5:
            taxa = list(range(min(n, 5)))
        genes.append(rand_gene(taxa, rng, p_dashed))

    xn = S.t2n()[0]                      # species node of taxon 0
    xtax = 0
    genes = [g for g in genes if xtax in g.t2n()]   # gene may lack the taxon
    assert genes
    Sm, _ = detach(S, xn)

    internal = [m for m in range(len(Sm.adj))
                if not Sm.is_leaf(m) and len(Sm.adj[m]) == 3]
    dirinfo = {}
    for m in internal:
        dirs = all_dirs(Sm, m)
        assert len(dirs) == 3
        dirinfo[m] = dirs

    # -------- L1 vs L2
    for m in internal:
        Ds = [dirinfo[m][d][1] for d in range(3)]
        for d in range(3):
            v1 = sum(T_direct(g, g.t2n()[xtax], Ds, d) for g in genes)
            v2 = sum(T_formula(g, g.t2n()[xtax], Ds, d) for g in genes)
            if v1 != v2:
                return f"FAIL L1!=L2 seed={seed} m={m} d={d}: direct={v1} formula={v2}"

    # -------- L4
    for g in genes:
        r = check_uniformity(g, g.t2n()[xtax], None)
        if r:
            return f"FAIL L4 seed={seed}: {r}"

    # -------- L5: one sweep per (gene, m, d) reproduces T for EVERY leaf
    for g in genes:
        for m in internal:
            Ds = [dirinfo[m][d][1] for d in range(3)]
            for d in range(3):
                seen = set()
                for ytax in sorted(g.t2n()):
                    yn = g.t2n()[ytax]
                    gamma = solid_cluster(g, yn)
                    if gamma in seen:
                        continue
                    seen.add(gamma)
                    got = sweep_T(g, Ds, d, gamma, 0)
                    for t, v in got.items():
                        want = T_formula(g, g.t2n()[t], Ds, d)
                        if v != want:
                            return (f"FAIL L5 seed={seed} m={m} d={d} "
                                    f"tax={t}: sweep={v} formula={want}")

    # -------- L3
    T = {}
    for m in internal:
        Ds = [dirinfo[m][d][1] for d in range(3)]
        for d in range(3):
            T[(m, d)] = sum(T_formula(g, g.t2n()[xtax], Ds, d) for g in genes)

    def G(R):
        u, v = R
        tot = 0
        for m in internal:
            par = None
            for idx, (w, tx) in enumerate(dirinfo[m]):
                if u in comp_nodes(Sm, m, w) or v in comp_nodes(Sm, m, w):
                    par = idx
                    break
            assert par is not None
            tot += T[(m, par)]
        return tot

    edges = [(i, w) for i in range(len(Sm.adj)) for w, _ in Sm.adj[i] if i < w]
    pi0 = edges[0]
    b0 = score_brute(attach(Sm, xtax, *pi0), genes)
    G0 = G(pi0)
    for pi in edges:
        b = score_brute(attach(Sm, xtax, *pi), genes)
        f = G(pi)
        if (b - b0) != (f - G0):
            return (f"FAIL L3 seed={seed} pi={pi}: brute_delta={b-b0} "
                    f"reroot_delta={f-G0}")
    return None


def main():
    nrun = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    nfail = 0
    for seed in range(nrun):
        r = check_instance(seed)
        if r:
            print(r)
            nfail += 1
    print(f"=== {nrun - nfail}/{nrun} instances passed ===")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
