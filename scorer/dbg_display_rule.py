#!/usr/bin/env python3
"""显示判据自检：SPEC 当前判据 vs 逐边扫描真值（随机多叉/S-D 基因树）。

SPEC 判据（u=LCA(a,b), v=LCA(c,d)）：
    ab|cd 显示 ⟺ (c,d ∉ subtree(u)) 或 (a,b ∉ subtree(v))
"""
import itertools
import random
import sys

from ref_brute import parse_newick, leaves_under, GeneTree
from difftest import rand_tree


def rule_displayed(g, a, b, c, d):
    u, v = g.lca(a, b), g.lca(c, d)
    su = set(leaves_under(g.ch, g.nm, u))
    sv = set(leaves_under(g.ch, g.nm, v))
    return (c not in su and d not in su) or (a not in sv and b not in sv)


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    rng = random.Random(seed)
    bad = checked = 0
    for _ in range(200):
        n = rng.randint(4, 9)
        taxa = [str(i) for i in range(n)]
        g = GeneTree(rand_tree(taxa, rng.choice([2, 3, 4]), rng.choice([0.0, 0.3, 0.5]), rng))
        leaves = sorted(g.taxa)
        for quad in itertools.combinations(leaves, 4):
            for p1, p2 in (((quad[0], quad[1]), (quad[2], quad[3])),
                           ((quad[0], quad[2]), (quad[1], quad[3])),
                           ((quad[0], quad[3]), (quad[1], quad[2]))):
                checked += 1
                if g.displayed(p1[0], p1[1], p2[0], p2[1]) != rule_displayed(g, p1[0], p1[1], p2[0], p2[1]):
                    bad += 1
                    if bad <= 3:
                        print("MISMATCH", p1, p2, "u,v =", g.lca(*p1), g.lca(*p2))
    print(f"{checked} 个四叶树检查完毕，与逐边真值不一致 {bad} 个")
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
