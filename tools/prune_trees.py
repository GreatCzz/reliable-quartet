#!/usr/bin/env python3
"""按叶子集合裁剪 Newick 树（物种树或基因树）。

用法:
  python3 prune_trees.py taxa N SEED OUT TREE      # 从 TREE 里随机选 N 个类群写到 OUT
  python3 prune_trees.py prune KEEP IN OUT         # 把 IN 里的每棵树裁到 KEEP 的叶子集合

基因树裁剪会压掉单叉节点，并把被压路径上的 D 标号合并到新边上。
"""
import sys, random

def parse(s):
    s = s.strip().rstrip(';')
    ch, lab = [], []
    leafLab = {}
    def rec(i):
        u = len(ch); ch.append([]); lab.append(''); i += 1
        while True:
            if s[i] == '(':
                c, i = rec(i); ch[u].append(c)
            else:
                j = i
                while j < len(s) and s[j] not in '(),': j += 1
                tok = s[i:j]
                # 叶子也可能带 S/D 标号（如 6670D）：拆出名字，标号记在叶子的 lab 上
                name, elab = tok, ''
                if len(tok) > 1 and tok[-1] in 'SD' and tok[-1].isdigit() is False:
                    name, elab = tok[:-1], tok[-1]
                leaf = len(ch); ch.append([]); lab.append(name); ch[u].append(leaf)
                leafLab[leaf] = elab
                i = j
            if s[i] == ',': i += 1; continue
            if s[i] == ')': i += 1; break
        j = i
        while j < len(s) and s[j] not in '(),;': j += 1
        lab[u] = s[i:j]; i = j
        return u, i
    root, _ = rec(0)
    return ch, lab, root

def prune_to(ch, lab, root, keep):
    # rec 返回 (子树 newick, 该子树根边是否 solid, 该子树是否已退化成单个叶子)
    def rec(u):
        if not ch[u]:
            return (lab[u], True, True) if lab[u] in keep else None
        kids = []
        for c in ch[u]:
            r = rec(c)
            if r is None: continue
            sub, solid, _ = r
            solid = solid and (lab[c] != 'D')
            kids.append((sub, solid, r[2]))
        if not kids: return None
        if len(kids) == 1:                     # 单叉节点被抑制，原样把孩子的结果顶上去
            return kids[0]
        # 叶子边一律 solid（题目口径：四条末端枝全实线）。绝不能把标号拼到叶名后面：
        # "1129" + "D" 会变成类群 "1129D"，在别的工具里就凭空多出一个假类群。
        parts = [sub + ('' if (solid or isLeaf) else 'D') for sub, solid, isLeaf in kids]
        return ('(' + ','.join(parts) + ')', True, False)
    r = rec(root)
    return None if r is None else r[0] + ';'

def main():
    mode = sys.argv[1]
    if mode == 'taxa':
        n, seed, out, src = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
        names = set()
        for line in open(src):
            if not line.strip(): continue
            ch, lab, root = parse(line)
            names = {lab[u] for u in range(len(ch)) if not ch[u]}
            break
        rng = random.Random(seed)
        S = rng.sample(sorted(names, key=int), n) if n < len(names) else sorted(names)
        open(out, 'w').write('\n'.join(S) + '\n')
        print(f'picked {len(S)} taxa', file=sys.stderr)
        return
    keep = {x.strip() for x in open(sys.argv[2]) if x.strip()}
    n = 0
    with open(sys.argv[3]) as f, open(sys.argv[4], 'w') as g:
        for line in f:
            line = line.strip()
            if not line: continue
            ch, lab, root = parse(line)
            r = prune_to(ch, lab, root, keep)
            if r: g.write(r + '\n'); n += 1
    print(f'pruned {n} trees -> {sys.argv[4]}', file=sys.stderr)

if __name__ == '__main__':
    main()
