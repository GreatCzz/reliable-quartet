#!/usr/bin/env python3
"""小工具：把树限制到给定类群集合（删叶 + 抑制单子节点，D 主导合并边标号）。"""


def restrict_tree(children, label, name, root, keep):
    new_ch, new_lab, new_nm = [], [], []

    def new_node(nm=None, lab=None):
        new_ch.append([])
        new_lab.append(lab)
        new_nm.append(nm)
        return len(new_ch) - 1

    def build(u):
        if not children[u]:
            if name[u] not in keep:
                return -1
            return new_node(name[u], label[u])
        kids = [build(c) for c in children[u]]
        kids = [k for k in kids if k >= 0]
        if not kids:
            return -1
        if len(kids) == 1:
            child = kids[0]
            lab = 'D' if (label[u] == 'D' or new_lab[child] == 'D') else (label[u] or new_lab[child])
            new_lab[child] = lab
            return child
        node = new_node(None, label[u])
        new_ch[node] = kids
        return node

    r = build(root)
    return new_ch, new_lab, new_nm, r


def render(children, label, name, u):
    if not children[u]:
        return name[u]
    return '(' + ','.join(render(children, label, name, c) for c in children[u]) + ')' + (label[u] or '')
