"""Well-formed adversarial instance generator.

Rules enforced (these are exactly what burned us):
  * every internal node has >= 1 child   (no empty "()" nodes)
  * every internal edge carries an EXPLICIT S or D label (real data does too)
  * no taxon appears twice in one tree
  * output validated by re-parsing with a strict checker
"""
import random, sys, subprocess, re, os

ADV = "/home/great_czz/project/演化算法助教工作/算法评分赛"

class Tree:
    def __init__(s): s.ch=[]; s.lab=[]; s.tax=[]
    def new(s,t=None,l=''):
        s.ch.append([]); s.lab.append(l); s.tax.append(t); return len(s.ch)-1
    def leaves(s): return [u for u in range(len(s.ch)) if s.tax[u] is not None]

def rand_binary(n,rng):
    """unrooted binary species tree, taxa 0..n-1"""
    t=Tree(); roots=[t.new(i) for i in range(n)]
    rng.shuffle(roots)
    while len(roots)>1:
        a=roots.pop(rng.randrange(len(roots))); b=roots.pop(rng.randrange(len(roots)))
        u=t.new(); t.ch[u]=[a,b]; roots.append(u)
    r=roots[0]
    # if the last node has 2 children it is a normal binary node; that's fine as newick
    return t, r

def rand_gene(taxa,rng,p_dashed,star=0,unary=False):
    t=Tree(); roots=[t.new(a) for a in taxa]
    rng.shuffle(roots)
    if star and len(roots)>=star:
        grp=[roots.pop(rng.randrange(len(roots))) for _ in range(star)]
        u=t.new(); t.ch[u]=grp; roots.append(u)
    while len(roots)>1:
        k=min(len(roots), rng.choice([2,2,2,3,4]))
        grp=[roots.pop(rng.randrange(len(roots))) for _ in range(k)]
        u=t.new(); t.ch[u]=grp; roots.append(u)
    r=roots[0]
    if unary:                      # a genuine single-child node, like rep01 has
        u=t.new(); t.ch[u]=[r]; r=u
    for u in range(len(t.ch)):
        if t.ch[u]:                # internal edge -> EXPLICIT label
            t.lab[u] = 'S' if rng.random()>p_dashed else 'D'
    t.lab[r]=''                    # the root has no parent edge
    return t,r

def to_newick(t,r):
    out=[]
    def rec(u):
        if t.tax[u] is not None: out.append(str(t.tax[u])); return
        assert len(t.ch[u])>=1, "empty node"
        out.append('(')
        for i,c in enumerate(t.ch[u]):
            if i: out.append(',')
            rec(c)
            if t.tax[c] is None and t.lab[c]: out.append(t.lab[c])
        out.append(')')
    rec(r); out.append(';')
    return ''.join(out)

def check(s):
    assert '()' not in s and '(,' not in s and ',)' not in s, "malformed: empty child"
    d=0
    for c in s:
        if c=='(': d+=1
        elif c==')':
            d-=1
            assert d>=0
    assert d==0, "unbalanced"

def build(name,n=9,ngenes=3,p_dashed=0.3,p_missing=0.1,star=0,unary=False,seed=0,caterpillar=False):
    rng=random.Random(seed)
    if caterpillar:
        t=Tree(); roots=[t.new(i) for i in range(n)]
        r=t.new(); t.ch[r]=[roots[0],roots[1]]
        for x in roots[2:]:
            u=t.new(); t.ch[u]=[r,x]; r=u
        S,sr=t,r
    else:
        S,sr=rand_binary(n,rng)
    spnw=to_newick(S,sr); check(spnw)
    genes=[]
    for _ in range(ngenes):
        taxa=[a for a in range(n) if rng.random()>p_missing]
        if len(taxa)<4: taxa=list(range(4))
        g,gr=rand_gene(sorted(taxa),rng,p_dashed,star,unary)
        gs=to_newick(g,gr); check(gs)
        genes.append(gs)
    d=f"{ADV}/rscore/tools/adv_cases"
    os.makedirs(d,exist_ok=True)
    open(f"{d}/{name}.sp.nwk","w").write(spnw+"\n")
    open(f"{d}/{name}.g.trees","w").write("\n".join(genes)+"\n")
    return spnw, genes

def validate(name):
    d=f"{ADV}/rscore/tools/adv_cases"
    out=subprocess.run([f"{ADV}/rscore/solver/score_fast","-t",f"{d}/{name}.sp.nwk",
                        "-g",f"{d}/{name}.g.trees","-j","1"],capture_output=True,text=True)
    txt=out.stdout+out.stderr
    m=re.search(r'"n_unused_taxa":(\d+),"n_ignored_taxa":(\d+)',out.stdout)
    warn = ("absent" in txt) or ("never occur" in txt)
    return (m.group(1),m.group(2)) if m else None, warn

CASES=[
 ("c1_dominant",   dict(p_dashed=0.10)),
 ("c2_all_small",  dict(p_dashed=0.55)),
 ("c3_half",       dict(p_dashed=0.30)),
 ("c4_star",       dict(p_dashed=0.20, star=6)),
 ("c5_unary",      dict(p_dashed=0.25, unary=True)),
 ("c6_caterpillar",dict(p_dashed=0.25, caterpillar=True)),
 ("c7_missing",    dict(p_dashed=0.30, p_missing=0.35)),
 ("c8_nodash",     dict(p_dashed=0.0)),
]
for i,(nm,kw) in enumerate(CASES):
    build(nm,seed=200+i,**kw)
ok=True
for nm,_ in CASES:
    r,warn=validate(nm)
    tag = "OK" if (r and r==('0','0') and not warn) else "*** BAD ***"
    if tag!="OK": ok=False
    print(f"{nm:16s} n_unused/n_ignored={r} warn={warn}  {tag}")
print("generator clean" if ok else "generator still produces bad input")
