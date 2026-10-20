#!/usr/bin/env python3
"""在 rep01 的小规模子集上跑一遍暴力打分，检查语义与量级。"""
import gzip, random, itertools, sys
import quartet_score as qs
from quartet_score import GeneTree, score
from test_quartet_score import restrict, render

random.seed(0)
fname = "../rep01.trees.gz"
ntax, ntree = 24, 30

with gzip.open(fname, 'rt') as f:
    lines = [l.strip() for l in f if l.strip()][:ntree]

# 选一个在这 30 棵树里都出现的类群集合
sets = []
for line in lines:
    t = qs.parse_newick(line)
    sets.append({u.name for u in [t] + []} if False else set())
# 直接扫叶子标签
import re
pat = re.compile(r'(?<=[(,])(\d+)(?=[,)])')
sets = [set(pat.findall(l)) for l in lines]
common = set.intersection(*sets)
print("common taxa across", ntree, "trees:", len(common))
taxa = sorted(random.sample(sorted(common), ntax))
print("subset:", taxa[:8], "...")

# 用一棵随机二叉物种树作为候选 S
def random_binary(labels):
    nodes = [[x] for x in labels]
    random.shuffle(nodes)
    while len(nodes) > 1:
        i = random.randrange(len(nodes)); a = nodes.pop(i)
        j = random.randrange(len(nodes)); b = nodes.pop(j)
        nodes.append(a + b)
    def rec(group):
        if len(group) == 1: return group[0]
        k = len(group) // 2
        return "(" + rec(group[:k]) + "," + rec(group[k:]) + ")"
    return rec(nodes[0]) + ";"

S = random_binary(taxa)
# 把基因树限制到这个子集
genes = []
for line in lines:
    t = qs.parse_newick(line)
    genes.append(render(restrict(t, set(taxa))) + ";")

import time
t0 = time.time()
total, per = score(S, genes)
print("score (reliable quartets) =", total, "per-gene first 5:", per[:5], "time %.1fs" % (time.time()-t0))
print("quartets in S:", len(qs.species_quartets(S)))
