#!/usr/bin/env python3
"""逐组合诊断：把"公式的显示判据"与暴力显示判据逐个四叶树比对。

对每个可靠的 (a,b|c,d) 组合，计算：
  disp_true  = 暴力边分裂判据
  disp_rule  = 公式的三类排除规则（锚点重合 / 祖先 + 叶子落在通往对方的孩子里）
若两者不一致，打印明细。这样能把"公式错"和"实现错"分开。
"""
import itertools
import random
import sys

from ref_brute import parse_newick, leaves_under, GeneTree


def is_anc(par, x, y):
    """y 是否为 x 的真祖先。"""
    cur = par[x]
    while cur != -1:
        if cur == y:
            return True
        cur = par[cur]
    return False


def child_toward(par, anc, desc):
    """anc 的孩子中，位于通向 desc 路径上的那个。"""
    cur = desc
    while par[cur] != anc:
        cur = par[cur]
    return cur


def rule_displayed(g, a, b, c, d):
    """按公式推导出的显示判据。"""
    u = g.lca(a, b)
    v = g.lca(c, d)
    if u == v:
        return False
    if is_anc(g.par, u, v):                     # v 是 u 的祖先?
        ch = child_toward(g.par, v, u)          # v 的、通向 u 的孩子
        if a in set(leaves_under(g.ch, g.nm, ch)) or b in set(leaves_under(g.ch, g.nm, ch)):
            return False
        return True
    if is_anc(g.par, v, u):                     # u 是 v 的祖先?
        ch = child_toward(g.par, u, v)
        if c in set(leaves_under(g.ch, g.nm, ch)) or d in set(leaves_under(g.ch, g.nm, ch)):
            return False
        return True
    return True                                  # 不可比 → 显示


def check_case(species_newick, gene_newick, verbose=True):
    sch, _, snm, _, _ = parse_newick(species_newick)
    g = GeneTree(gene_newick)
    bad = 0
    for v in range(len(sch)):
        if not sch[v] or len(sch[v]) != 2:
            continue
        L = sorted(set(leaves_under(sch, snm, sch[v][0])) & g.taxa)
        R = sorted(set(leaves_under(sch, snm, sch[v][1])) & g.taxa)
        for a, b in itertools.combinations(L, 2):
            for c, d in itertools.combinations(R, 2):
                u, v2 = g.lca(a, b), g.lca(c, d)
                solid = (g.path_solid(a, u) and g.path_solid(b, u)
                         and g.path_solid(c, v2) and g.path_solid(d, v2))
                if not solid:
                    continue
                t = g.displayed(a, b, c, d)
                r = rule_displayed(g, a, b, c, d)
                if t != r:
                    bad += 1
                    if verbose:
                        print(f"  S节点 {v}: {a},{b}|{c},{d} u={u} v={v2} "
                              f"displayed_true={t} rule={r}")
    return bad


def rand_tree(taxa, maxdeg=2, dprob=0.3, rng=None):
    rng = rng or random
    def rec(items):
        if len(items) == 1:
            return str(items[0])
        k = rng.choice([2] * 4 + list(range(3, maxdeg + 1))) if maxdeg > 2 else 2
        k = max(2, min(k, len(items)))
        items = items[:]
        rng.shuffle(items)
        parts = [items[i::k] for i in range(k)]
        parts = [p for p in parts if p]
        if len(parts) == 1:
            return str(parts[0][0])
        return '(' + ','.join(rec(p) for p in parts) + ')' + ('D' if rng.random() < dprob else 'S')
    return rec(list(taxa)) + ';'


if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    tested = 0
    for trial in range(500):
        n = rng.randint(4, 9)
        taxa = [str(i) for i in range(n)]
        S = rand_tree(taxa, 2, 0.0, rng)
        g = rand_tree(taxa, rng.choice([2, 3, 4]), 0.35, rng)
        bad = check_case(S, g, verbose=(tested < 5))
        if bad:
            tested += 1
            print(f"[case {trial}] S={S}\n           g={g}  不一致数量={bad}\n")
            if tested >= 5:
                break
    if tested == 0:
        print("500 组随机用例：显示判据与暴力完全一致 → 公式的显示规则没问题，问题在计数实现")
