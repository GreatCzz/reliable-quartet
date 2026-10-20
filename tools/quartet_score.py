#!/usr/bin/env python3
"""可靠四叶树打分原型（暴力版，仅用于小规模验证语义）。

定义（依据 week5.pdf p44-47 的可靠三叶树定义推广）：
- 基因树 g 为带 S/D 边标的（可非二叉）树；叶子边无标号，视为实线。
- 四叶树 ab|cd 在 g 中"显示"：存在一条内部边 e=(p,x)，使得 a,b 在 x 侧、
  c,d 在 p 侧（或反之）。
- ab|cd 在 g 中"可靠"：存在这样的 e，且 p->a、p->b、x->c、x->d 四条末端枝
  上的每条边都是实线（S 或无标号）；锚点之间的那条边 e 本身可以是虚边 D。
- 物种树 S 与 g 的公共可靠四叶树 = S 显示的且 g 中可靠的四叶树个数。
"""
import itertools
import sys

class Node:
    __slots__ = ("name", "label", "children", "parent")
    def __init__(self, name=None, label=None):
        self.name = name; self.label = label; self.children = []; self.parent = None

def parse_newick(s):
    """Parse newick with optional S/D label after ')' -> edge (parent->node) label."""
    s = s.strip().rstrip(';')
    i, n = 0, len(s)
    root = None
    stack = []
    while i < n:
        c = s[i]
        if c == '(':
            node = Node()
            if stack: node.parent = stack[-1]; stack[-1].children.append(node)
            else: root = node
            stack.append(node); i += 1
        elif c == ')':
            j = i + 1; lab = ''
            while j < n and s[j] not in '(),':
                lab += s[j]; j += 1
            if stack:
                stack[-1].label = lab or None
            i = j; stack.pop()
        elif c == ',':
            i += 1
        else:
            j = i; tok = ''
            while j < n and s[j] not in '(),':
                tok += s[j]; j += 1
            leaf = Node(tok)
            leaf.parent = stack[-1]; stack[-1].children.append(leaf)
            i = j
    return root

def leaves_of(root):
    out = {}
    def rec(u):
        if not u.children: out[u] = {u.name}; return out[u]
        s = set()
        for c in u.children:
            s |= rec(c)
        out[u] = s
        return s
    rec(root)
    return out

def is_solid(u):
    """edge above u is solid?"""
    return u.label in (None, 'S')

def path_all_solid(u, ancestor):
    """edges on path from ancestor (exclusive) down to u (inclusive) all solid?"""
    cur = u
    while cur is not ancestor:
        if not is_solid(cur): return False
        cur = cur.parent
    return True

class GeneTree:
    def __init__(self, newick):
        self.root = parse_newick(newick)
        self.leafsets = leaves_of(self.root)
        self.nodes = list(self.leafsets)
        self.taxa = {u.name for u in self.nodes if not u.children}
        # internal edges (p, x)
        self.edges = [(u.parent, u) for u in self.nodes if u.children and u.parent is not None]

    def quartet_state(self, a, b, c, d):
        """Return (displayed, reliable) for quartet ab|cd.

        displayed: some internal edge e splits {a,b} from {c,d}
        reliable: displayed, and the four pendant paths from u=LCA(a,b) and
        v=LCA(c,d) down to a,b,c,d consist only of solid edges.
        """
        side = {a, b}; other = {c, d}
        displayed = False
        for p, x in self.edges:
            ls = self.leafsets[x]
            if (side <= ls and not (other & ls)) or (other <= ls and not (side & ls)):
                displayed = True
                break
        if not displayed:
            return False, False
        u = self.lca(a, b); v = self.lca(c, d)
        if u is v:
            return True, False
        ok = (path_all_solid(self._leaf(a), u) and path_all_solid(self._leaf(b), u)
              and path_all_solid(self._leaf(c), v) and path_all_solid(self._leaf(d), v))
        return True, ok

    def lca(self, x, y):
        nx, ny = self._leaf(x), self._leaf(y)
        seen = set()
        cur = nx
        while cur is not None:
            seen.add(id(cur)); cur = cur.parent
        cur = ny
        while id(cur) not in seen:
            cur = cur.parent
        return cur

    def _leaf(self, name):
        for u in self.nodes:
            if u.name == name and not u.children:
                return u
        raise KeyError(name)

def species_quartets(S_newick):
    """Return dict quartet-frozenset-ab|cd -> None; S must be binary."""
    root = parse_newick(S_newick)
    ls = leaves_of(root)
    edges = [(u.parent, u) for u in ls if u.children and u.parent is not None]
    quartets = set()
    # every quartet of S is displayed by exactly one edge in a binary tree:
    # the one separating the two pairs
    all_leaves = {u.name for u in ls if not u.children}
    for p, x in edges:
        A = ls[x]; B = all_leaves - A
        if len(A) < 2 or len(B) < 2: continue
        for a, b in itertools.combinations(sorted(A), 2):
            for c, d in itertools.combinations(sorted(B), 2):
                quartets.add(frozenset((frozenset((a, b)), frozenset((c, d)))))
    return quartets

def score(S_newick, gene_newicks, verbose=False):
    S = species_quartets(S_newick)
    total = 0
    per_gene = []
    for gn in gene_newicks:
        g = GeneTree(gn)
        cnt = 0
        for q in S:
            (p1, p2) = tuple(q)
            a, b = sorted(p1); c, d = sorted(p2)
            if not set((a, b, c, d)) <= g.taxa: continue
            disp, rel = g.quartet_state(a, b, c, d)
            if rel: cnt += 1
        per_gene.append(cnt)
        total += cnt
    if verbose: print("per-gene:", per_gene[:10], "...")
    return total, per_gene

if __name__ == "__main__":
    g1 = "((a,b)S,(c,d)S);"
    g2 = "((a,b)S,(c,d)D);"
    g3 = "((a,b)S,c,d);"           # polytomy: ab|cd not displayed
    g4 = "(((a,b)S,c)S,d);"        # ab|cd displayed via internal edges
    S  = "((a,b),(c,d));"
    for name, g in [("g1", g1), ("g2", g2), ("g3", g3), ("g4", g4)]:
        t, per = score(S, [g])
        print(f"{name}: reliable quartets with S = {t}")
