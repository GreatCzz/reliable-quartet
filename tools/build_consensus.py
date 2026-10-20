#!/usr/bin/env python3
"""基线求解器 v0：把基因树中反复出现的子簇汇总成共识层级，再二值化成物种树。

依据数据画像：reliable 四叶树几乎全部来自已解析的小簇（樱桃/三叶簇/…），
因此先搭"簇骨架"，再把剩余叶子二值化挂上去。

用法: python3 build_consensus.py ../rep01.trees.gz out.tree [min_count]
"""
import gzip, sys, collections, random
sys.setrecursionlimit(200000)

def parse_fast(s):
    names = []; children = []; labels = []
    def nn(n=None, lab=None):
        names.append(n); children.append([]); labels.append(lab); return len(names) - 1
    st = [[]]; i = 0; n = len(s)
    while i < n:
        c = s[i]
        if c == '(':
            st.append([]); i += 1
        elif c == ')':
            kids = st.pop(); j = i + 1; lab = ''
            while j < n and s[j] not in '(),;':
                lab += s[j]; j += 1
            k = nn(None, lab); children[k] = kids; st[-1].append(k); i = j
        elif c in ',;':
            i += 1
        else:
            j = i; t = ''
            while j < n and s[j] not in '(),;':
                t += s[j]; j += 1
            st[-1].append(nn(t)); i = j
    return names, children, labels

def clades_of(s, max_size=5000):
    """返回 [(clade_tuple, subtree_all_solid)]，不含单叶与全树。"""
    names, children, labels = parse_fast(s)
    n = len(names)
    childset = set()
    for x in children:
        childset.update(x)
    root = next(i for i in range(n) if i not in childset)
    sets = [None] * n
    solid = [True] * n
    out = []
    for u in range(n):                      # 子节点 id 一定小于父节点
        if not children[u]:
            sets[u] = (names[u],); solid[u] = True
        else:
            parts = []
            ok = labels[u] in (None, 'S')
            for c in children[u]:
                parts.extend(sets[c]); ok = ok and solid[c]
            sets[u] = tuple(parts); solid[u] = ok
            if u != root and 2 <= len(sets[u]) <= max_size:
                out.append((tuple(sorted(sets[u])), ok))
    taxa = {names[i] for i in range(n) if not children[i]}
    return out, taxa

class N:
    __slots__ = ('kids', 'name', '_set')
    def __init__(self, name=None, kids=None):
        self.kids = kids if kids is not None else []
        self.name = name
        self._set = None
    def leafset(self):
        if self.name is not None:
            return frozenset((self.name,))
        if self._set is None:
            s = set()
            for k in self.kids:
                s |= k.leafset()
            self._set = frozenset(s)
        return self._set
    def invalidate(self):
        self._set = None

def build_tree(clade_counts, all_taxa, min_count):
    root = N(kids=[N(name=t) for t in sorted(all_taxa)])
    items = sorted(clade_counts.items(), key=lambda kv: (-kv[1], -len(kv[0])))
    inserted = 0
    for clade, cnt in items:
        if cnt < min_count:
            continue
        C = frozenset(clade)
        cur = root
        while True:
            nxt = None
            for k in cur.kids:
                if C <= k.leafset():
                    nxt = k; break
            if nxt is None or nxt.name is not None:
                break
            cur = nxt
        if cur.name is not None:
            continue
        inside = [k for k in cur.kids if k.leafset() <= C]
        outside = [k for k in cur.kids if not (k.leafset() & C)]
        if len(inside) < 2 or len(inside) + len(outside) != len(cur.kids):
            continue                      # 冲突 / 已存在 / 跨越多个孩子
        newnode = N(kids=inside)
        cur.kids = outside + [newnode]
        cur.invalidate()
        inserted += 1
    return root, inserted

def binarize(node, rng):
    """平衡二值化：每轮把相邻两个孩子合并，深度 O(log k)。"""
    if node.name is not None:
        return node
    kids = [binarize(k, rng) for k in node.kids]
    rng.shuffle(kids)
    while len(kids) > 1:
        nxt = []
        for i in range(0, len(kids) - 1, 2):
            nxt.append(N(kids=[kids[i], kids[i + 1]]))
        if len(kids) % 2:
            nxt.append(kids[-1])
        kids = nxt
    return kids[0]

def render(node):
    if node.name is not None:
        return node.name
    return '(' + ','.join(render(k) for k in node.kids) + ')'

def main():
    fname = sys.argv[1] if len(sys.argv) > 1 else "../rep01.trees.gz"
    out = sys.argv[2] if len(sys.argv) > 2 else "/tmp/team01_consensus.tree"
    min_count = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    counts = collections.Counter(); all_taxa = set()
    ntrees = 0
    with gzip.open(fname, 'rt') if fname.endswith('.gz') else open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            ntrees += 1
            cl, taxa = clades_of(line)
            all_taxa |= taxa
            for clade, _solid in cl:
                counts[clade] += 1
    print(f"trees: {ntrees}, taxa: {len(all_taxa)}, distinct clades: {len(counts)}")
    tree, inserted = build_tree(counts, all_taxa, min_count)
    print("clusters inserted:", inserted)
    tree = binarize(tree, random.Random(0))
    newick = render(tree) + ";"
    with open(out, "w") as f:
        f.write(newick + "\n")
    print("wrote", out, "leaves:", newick.count(',') + 1, "len", len(newick))

if __name__ == "__main__":
    main()
