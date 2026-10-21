#!/usr/bin/env bash
# 一键自检：构建 + 打分器对拍 + 提交树合法性 + C++/Python 一致性
set -uo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
fail=0
pass() { echo "  [PASS] $*"; }
bad()  { echo "  [FAIL] $*"; fail=1; }

echo "== 1/4 编译 solver =="
if ( cd "$ROOT/rscore/solver" && g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src solver.cpp -o solver ); then
  pass "solver 编译通过"
else
  bad "solver 编译失败"; exit 1
fi

echo "== 2/4 打分器对拍（Python 暴力 vs 公式） =="
if ( cd "$ROOT/scorer" && python3 difftest.py 200 0 >/tmp/verify_difftest.log 2>&1 ); then
  pass "$(tail -1 /tmp/verify_difftest.log)"
else
  bad "对拍失败（见 /tmp/verify_difftest.log）"
fi

echo "== 3/4 提交树合法性（叶子/标签/二叉性） =="
python3 - "$ROOT" <<'PYCHECK'
import glob, sys, os
root = sys.argv[1]

def parse(s):
    """返回 (每个内部节点的孩子数列表, 叶子标签列表)；结构非法则抛异常。"""
    s = s.strip().rstrip(';')
    deg, labels = [], []
    stack = [[]]
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == '(':
            stack.append([]); i += 1
        elif c == ')':
            kids = stack.pop()
            if not kids:
                raise ValueError("empty internal node")
            deg.append(len(kids))
            stack[-1].append(None)
            i += 1
            while i < n and s[i] not in '(),;':
                i += 1
        elif c == ',':
            i += 1
        else:
            j = i
            while j < n and s[j] not in '(),;':
                j += 1
            lab = s[i:j].strip()
            if not lab:
                raise ValueError("empty leaf label")
            labels.append(lab)
            stack[-1].append(lab)
            i = j
    if len(stack) != 1 or len(stack[0]) != 1:
        raise ValueError("malformed newick")
    return deg, labels

ok = True
for f in sorted(glob.glob(os.path.join(root, 'submissions', '*.tree'))):
    try:
        deg, labels = parse(open(f).read())
        binary = all(d == 2 for d in deg)
        unique = len(set(labels)) == len(labels)
        full = (len(labels) == 10000 and unique
                and min(map(int, labels)) == 0 and max(map(int, labels)) == 9999)
        good = binary and full
        msg = f"leaves={len(labels)} binary={binary} unique={unique}"
    except Exception as e:
        good, msg = False, f"parse error: {e}"
    print(f"  [{'PASS' if good else 'FAIL'}] {os.path.basename(f)}: {msg}")
    ok = ok and good
sys.exit(0 if ok else 1)
PYCHECK
[ $? -eq 0 ] || fail=1

echo "== 4/4 C++ 与 Python 打分一致性（40 类群真实数据） =="
if [ -f "$ROOT/rep01.trees.gz" ]; then
  tmp=$(mktemp -d)
  python3 - "$ROOT" "$tmp" <<'PY'
import gzip, random, re, sys, os
root, tmp = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(root, 'scorer'))
from ref_brute import parse_newick
from test_utils import restrict_tree, render
with gzip.open(os.path.join(root, 'rep01.trees.gz'), 'rt') as f:
    lines = [l.strip() for l in f if l.strip()][:30]
pat = re.compile(r'(?<=[(,])(\d+)(?=[,)])')
common = set.intersection(*[set(pat.findall(l)) for l in lines])
rng = random.Random(0)
taxa = sorted(rng.sample(sorted(common), 40))
out = []
for l in lines:
    ch, lab, nm, par, rt = parse_newick(l)
    nc, nl, nn, nr = restrict_tree(ch, lab, nm, rt, set(taxa))
    out.append(render(nc, nl, nn, nr) + ";")
open(os.path.join(tmp, 'g40.trees'), 'w').write("\n".join(out) + "\n")
def rnd(items):
    if len(items) == 1: return items[0]
    rng.shuffle(items); k = len(items)//2
    return "(" + rnd(items[:k]) + "," + rnd(items[k:]) + ")"
open(os.path.join(tmp, 's40.tree'), 'w').write(rnd(taxa[:]) + ";\n")
PY
  cpp=$("$ROOT/rscore/solver/solver" -s "$tmp/s40.tree" -i "$tmp/g40.trees" -t 4 --sweeps 0 2>&1 | grep "initial score" | sed 's/.*= //' | awk '{print $1}')
  py=$(cd "$ROOT/scorer" && python3 - "$tmp" <<'PY'
import sys
from fast_formula import fast_score
tmp = sys.argv[1]
S = open(tmp + '/s40.tree').read().strip()
genes = [l.strip() for l in open(tmp + '/g40.trees') if l.strip()]
print(sum(fast_score(S, g) for g in genes))
PY
)
  if [ "$cpp" = "$py" ]; then pass "C++=$cpp == Python=$py"; else bad "C++=$cpp != Python=$py"; fi
  rm -rf "$tmp"
else
  echo "  [SKIP] 没有 rep01.trees.gz"
fi

echo
if [ $fail -eq 0 ]; then echo "全部通过 ✅"; else echo "存在失败项 ❌"; fi
exit $fail
