#!/usr/bin/env python3
"""差分测试：快速公式版 vs 暴力参考版。

用法：python3 difftest.py [轮数] [随机种子]
"""
import random
import sys

from ref_brute import brute_score
from fast_formula import fast_score_detailed


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


def run(rounds=500, seed=0, nmin=4, nmax=9, verbose=True):
    rng = random.Random(seed)
    bad = 0
    for t in range(rounds):
        n = rng.randint(nmin, nmax)
        taxa = [str(i) for i in range(n)]
        S = rand_tree(taxa, 2, 0.0, rng)
        g = rand_tree(taxa, rng.choice([2, 3, 4]), rng.choice([0.0, 0.2, 0.35, 0.5]), rng)
        exp, exp_per = brute_score(S, g)
        got, got_per = fast_score_detailed(S, g)
        if got != exp:
            bad += 1
            if verbose and bad <= 3:
                print(f"[MISMATCH] S={S}\n           g={g}\n           期望 {exp}  实际 {got}")
                for v in exp_per:
                    if exp_per[v] != got_per.get(v, 0):
                        print(f"             节点 {v}: 期望 {exp_per[v]} 实际 {got_per.get(v, 0)}")
    print(f"{rounds} 组随机用例（n={nmin}..{nmax}）：不一致 {bad} 组")
    return bad


if __name__ == '__main__':
    rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    sys.exit(1 if run(rounds, seed) else 0)
