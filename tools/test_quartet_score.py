#!/usr/bin/env python3
"""随机小树交叉验证 quartet_score 的 displayed/reliable 判定。

independent 实现：先"删掉无关叶子 + 抑制单子节点"得到诱导 4 叶树，
再在诱导树上按候选配对判断分裂与末端枝实线情况。
"""
import random, itertools, sys
import quartet_score as qs
from quartet_score import GeneTree

def random_tree_str(taxa):
    def rec(leaves):
        if len(leaves) == 1:
            return str(leaves[0])
        k = min(random.choice([2, 2, 2, 3]), len(leaves))
        shuffled = leaves[:]
        random.shuffle(shuffled)
        parts = [shuffled[i::k] for i in range(k)]
        parts = [p for p in parts if p]
        if len(parts) == 1:
            return str(parts[0][0])
        lab = random.choice(['S', 'S', 'S', 'D'])
        return "(" + ",".join(rec(p) for p in parts) + ")" + lab
    return rec(list(taxa)) + ";"

def render(u):
    if not u.children:
        return u.name
    return "(" + ",".join(render(c) for c in u.children) + ")" + (u.label or "")

def restrict(root, keep):
    """重建式限制：只保留 keep 中的叶子，抑制单子节点，边标号 D 主导。"""
    def build(u):
        if not u.children:
            if u.name not in keep:
                return None
            return qs.Node(u.name, u.label)
        kids = [build(c) for c in u.children]
        kids = [k for k in kids if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            child = kids[0]
            lab = 'D' if (u.label == 'D' or child.label == 'D') else (u.label or child.label)
            child.label = lab
            child.parent = None
            return child
        node = qs.Node(None, u.label)
        for k in kids:
            k.parent = node
            node.children.append(k)
        return node
    return build(root)

def leaves_under(u):
    out, st = set(), [u]
    while st:
        v = st.pop()
        if v.children:
            st.extend(v.children)
        else:
            out.add(v.name)
    return out

def induced_expected(tree_str, p1, p2):
    """在诱导 4 叶树上判断：先看分裂是否为给定配对，再用 LCA 锚点查末端枝实线。"""
    root = restrict(qs.parse_newick(tree_str), set(p1) | set(p2))
    nodes, st = [], [root]
    while st:
        u = st.pop(); nodes.append(u); st.extend(u.children)
    leaves = [u for u in nodes if not u.children]
    if len(leaves) != 4 or not any(u.children for u in nodes):
        return (False, False)
    target1, target2 = set(p1), set(p2)
    allq = target1 | target2
    split_ok = False
    for x in nodes:
        p = x.parent
        if p is None:
            continue
        A = leaves_under(x)
        B = allq - A
        if (A == target1 and B == target2) or (A == target2 and B == target1):
            split_ok = True
            break
    if not split_ok:
        return (False, False)
    def leaf(name):
        return next(u for u in nodes if not u.children and u.name == name)
    def lca(x, y):
        seen, cur = set(), x
        while cur is not None:
            seen.add(id(cur)); cur = cur.parent
        cur = y
        while id(cur) not in seen:
            cur = cur.parent
        return cur
    def solid(u, anc):
        cur = u
        while cur is not anc:
            if cur.label == 'D':
                return False
            cur = cur.parent
        return True
    u1 = lca(leaf(p1[0]), leaf(p1[1]))
    u2 = lca(leaf(p2[0]), leaf(p2[1]))
    if u1 is u2:
        return (True, False)
    ok = all(solid(leaf(z), u1) for z in p1) and all(solid(leaf(z), u2) for z in p2)
    return (True, ok)

if __name__ == "__main__":
    random.seed(7)
    taxa = [str(i) for i in range(7)]
    mism = 0
    for trial in range(80):
        s = random_tree_str(taxa)
        g = GeneTree(s)
        for quad in itertools.combinations(taxa, 4):
            pairs = [((quad[0], quad[1]), (quad[2], quad[3])),
                     ((quad[0], quad[2]), (quad[1], quad[3])),
                     ((quad[0], quad[3]), (quad[1], quad[2]))]
            for p1, p2 in pairs:
                got = g.quartet_state(p1[0], p1[1], p2[0], p2[1])
                exp = induced_expected(s, p1, p2)
                if got != exp:
                    mism += 1
                    print("MISMATCH", s, p1, p2, "got", got, "exp", exp)
                    if mism > 5:
                        sys.exit(1)
    print("done, mismatches:", mism)
