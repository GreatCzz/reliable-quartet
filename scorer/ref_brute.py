#!/usr/bin/env python3
"""打分器参考实现（暴力，绝对可信）：公共可靠四叶树计数。

定义
----
- 基因树 g：可多叉；内部边带 S/D 标号（叶子边无标号，视为实线）。
- 四叶树 ab|cd 在 g 中"显示"：存在一条边（根以下的某个节点所对应的边），
  它的一侧恰好含 a,b 而不含 c,d（或反之）。
- ab|cd 在 g 中"可靠"：路径 a→b、c→d 上的每条边都是实线
  （等价于 LCA(a,b)、LCA(c,d) 到四个叶子的四条末端枝全实线）。
- 物种树 S 二叉：每个内部节点的两个孩子簇 l, r 贡献
  #{ab|cd : a,b ∈ l, c,d ∈ r, 且 q 在 g 中可靠}。
"""
import itertools


def parse_newick(s):
    """返回 (children, label, name, parent, root)；label 是"到父节点那条边"的标号。"""
    children, label, name = [], [], []

    def new_node(nm=None, lab=None):
        children.append([])
        label.append(lab)
        name.append(nm)
        return len(children) - 1

    stack = [[]]
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == '(':
            stack.append([])
            i += 1
        elif c == ')':
            kids = stack.pop()
            j = i + 1
            lab = ''
            while j < n and s[j] not in '(),;':
                lab += s[j]
                j += 1
            node = new_node(None, lab or None)
            children[node] = kids
            stack[-1].append(node)
            i = j
        elif c in ',;':
            i += 1
        else:
            j = i
            while j < n and s[j] not in '(),;':
                j += 1
            stack[-1].append(new_node(s[i:j]))
            i = j
    parent = [-1] * len(children)
    for u in range(len(children)):
        for c in children[u]:
            parent[c] = u
    child_set = {c for kids in children for c in kids}
    root = next(u for u in range(len(children)) if u not in child_set)
    return children, label, name, parent, root


def leaves_under(children, name, u):
    out, stack = [], [u]
    while stack:
        x = stack.pop()
        if not children[x]:
            out.append(name[x])
        else:
            stack.extend(children[x])
    return out


class GeneTree:
    def __init__(self, newick):
        self.ch, self.lab, self.nm, self.par, self.root = parse_newick(newick)
        self.leaf_node = {self.nm[u]: u for u in range(len(self.ch)) if not self.ch[u]}
        self.taxa = set(self.leaf_node)

    def solid_edge(self, u):
        """节点 u 到其父节点的边是否实线。"""
        return self.lab[u] != 'D'

    def path_solid(self, leaf_name, anc):
        cur = self.leaf_node[leaf_name]
        while cur != anc:
            if not self.solid_edge(cur):
                return False
            cur = self.par[cur]
        return True

    def lca(self, a, b):
        x, y = self.leaf_node[a], self.leaf_node[b]
        seen = set()
        while x != -1:
            seen.add(x)
            x = self.par[x]
        while y not in seen:
            y = self.par[y]
        return y

    def displayed(self, a, b, c, d):
        A, B = {a, b}, {c, d}
        for u in range(len(self.ch)):
            if u == self.root:
                continue
            ls = set(leaves_under(self.ch, self.nm, u))
            if (A <= ls and not (B & ls)) or (B <= ls and not (A & ls)):
                return True
        return False

    def reliable(self, a, b, c, d):
        if not self.displayed(a, b, c, d):
            return False
        u, v = self.lca(a, b), self.lca(c, d)
        return (self.path_solid(a, u) and self.path_solid(b, u)
                and self.path_solid(c, v) and self.path_solid(d, v))


def _edge_leafsets(sch, snm, sroot):
    return {u: set(leaves_under(sch, snm, u)) for u in range(len(sch))}


def brute_score(species_newick, gene_newick):
    """无根四叶树真值：每个四叶树只数一次（逐边分裂判定），与根位置无关。

    返回 (总分, {物种树节点: 该节点的空占位}) —— 第二个返回值仅为兼容旧接口。
    """
    sch, _, snm, _, sroot = parse_newick(species_newick)
    g = GeneTree(gene_newick)
    leafsets = _edge_leafsets(sch, snm, sroot)
    all_leaves = leafsets[sroot]
    taxa = sorted(set(all_leaves) & g.taxa)
    total = 0
    for quad in itertools.combinations(taxa, 4):
        for (p1, p2) in (((quad[0], quad[1]), (quad[2], quad[3])),
                         ((quad[0], quad[2]), (quad[1], quad[3])),
                         ((quad[0], quad[3]), (quad[1], quad[2]))):
            A, B = set(p1), set(p2)
            shown = False
            for u in range(len(sch)):
                if u == sroot:
                    continue
                Su = leafsets[u]
                Sv = all_leaves - Su
                if (A <= Su and B <= Sv) or (B <= Su and A <= Sv):
                    shown = True
                    break
            if shown and g.reliable(p1[0], p1[1], p2[0], p2[1]):
                total += 1
                break
    return total, {}
