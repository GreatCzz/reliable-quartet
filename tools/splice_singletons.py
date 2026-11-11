#!/usr/bin/env python3
"""去掉基因树里"只有一个孩子"的内部节点（splice singletons）。

ASTRAL-IV / ASTER2 的 huffman() 在单叉节点上会崩或 segfault，所以喂给它们之前必须 splice。

边标签合并规则与 `score_fast.cpp` 的 `splice1()` 一致：取优先级更高的那个
（D=3 > S/default=2 > C=1）。直觉：一条边只要有一段是虚边，整条就是虚边。

用法:
    python3 tools/splice_singletons.py <in.trees> <out.trees>

输入输出都是每行一棵 Newick；纯拓扑 + S/D 标签，无枝长。
"""
import sys


def parse(line):
    """返回 (children, labels, names, root)。children[i] 为孩子 id 列表。"""
    s = line.strip()
    children, labels, names = [], [], []
    st = []
    root = -1
    i, n = 0, len(s)

    def new(name="", lab=""):
        children.append([])
        labels.append(lab)
        names.append(name)
        return len(children) - 1

    while i < n:
        c = s[i]
        if c == '(':
            nd = new()
            if st:
                children[st[-1]].append(nd)
            else:
                root = nd
            st.append(nd)
            i += 1
        elif c == ',':
            i += 1
        elif c == ')':
            i += 1
            j = i
            while j < n and s[j] not in '(),;':
                j += 1
            if st:
                labels[st[-1]] = s[i:j]
                st.pop()
            i = j
        elif c == ';':
            i += 1
        else:
            j = i
            while j < n and s[j] not in '(),;':
                j += 1
            nd = new(name=s[i:j])
            if st:
                children[st[-1]].append(nd)
            else:
                root = nd
            i = j
    return children, labels, names, root


PRI = {'D': 3, 'S': 2, 'C': 1, '': 2}


def splice(children, labels, names, root):
    """返回新的 (children, labels, names, root)；只改 children 里的引用，节点数组原样保留。"""
    # 后序：先把孩子的单叉链压平
    order, st = [], [root]
    while st:
        t = st.pop()
        order.append(t)
        st.extend(children[t])
    order.reverse()

    changed = False
    for t in order:
        new_kids = []
        for c in children[t]:
            # 沿单叉链往下压，合并标签
            while len(children[c]) == 1:
                only = children[c][0]
                if PRI[labels[c]] >= PRI[labels[only]]:
                    labels[only] = labels[c]
                c = only
                changed = True
            new_kids.append(c)
        if new_kids != children[t]:      # 注意：必须比内容，不能比长度
            children[t] = new_kids
            changed = True

    # 根如果只有一个孩子，往下提一层
    while len(children[root]) == 1:
        only = children[root][0]
        if PRI[labels[root]] >= PRI[labels[only]]:
            labels[only] = labels[root]
        root = only
        changed = True
    return children, labels, names, root, changed


def render(children, labels, names, t):
    if not children[t]:
        return names[t]
    return "(" + ",".join(render(children, labels, names, c) for c in children[t]) + ")" + labels[t]


def main():
    if len(sys.argv) != 3:
        sys.stderr.write(__doc__)
        return 2
    src, dst = sys.argv[1], sys.argv[2]
    n_in = n_out = n_fix = 0
    with open(src) as fi, open(dst, 'w') as fo:
        for line in fi:
            line = line.strip()
            if not line:
                continue
            n_in += 1
            children, labels, names, root = parse(line)
            children, labels, names, root, changed = splice(children, labels, names, root)
            if changed:
                n_fix += 1
            fo.write(render(children, labels, names, root) + ";\n")
            n_out += 1
    sys.stderr.write(f"[splice] {n_in} 棵树读入，{n_out} 棵写出，{n_fix} 棵被改\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
