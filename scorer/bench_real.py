#!/usr/bin/env python3
"""在真实数据子集上测速并外推全量耗时。

用法: python3 bench_real.py [类群数] [基因树数]
"""
import gzip
import random
import re
import sys
import time

from fast_formula import fast_score_detailed
from ref_brute import parse_newick, leaves_under
from test_utils import restrict_tree, render

def load(n_taxa, n_trees, seed=0):
    fn = "../rep01.trees.gz"
    with gzip.open(fn, "rt") as f:
        lines = [ln.strip() for ln in f if ln.strip()][:n_trees]
    pat = re.compile(r'(?<=[(,])(\d+)(?=[,)])')
    common = set.intersection(*[set(pat.findall(l)) for l in lines])
    rng = random.Random(seed)
    taxa = sorted(rng.sample(sorted(common), n_taxa))
    ch, _, nm, _, _ = parse_newick(lines[0])
    from ref_brute import parse_newick as pn
    genes = []
    for line in lines:
        c, lab, name, par, root = pn(line)
        nc, nl, nn, nr = restrict_tree(c, lab, name, root, set(taxa))
        genes.append(render(nc, nl, nn, nr) + ";")
    return taxa, genes


if __name__ == "__main__":
    n_taxa = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    n_trees = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    sys.path.insert(0, ".")
    from difftest import rand_tree
    import random as _r
    taxa, genes = load(n_taxa, n_trees)
    rng = _r.Random(0)
    S = rand_tree(sorted(taxa), 2, 0.0, rng)
    t0 = time.time()
    total = 0
    for g in genes:
        total += fast_score_detailed(S, g)[0]
    dt = time.time() - t0
    print(f"类群 {n_taxa}、基因树 {n_trees}：总分 {total}，用时 {dt:.1f}s "
          f"（{dt / n_trees * 1000:.0f} ms/树）")
