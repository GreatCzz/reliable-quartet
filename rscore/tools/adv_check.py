"""对拍：sector_exact（C++, 真重算） vs Python 独立暴力参考（sector_math_check）"""
import sys, subprocess, re
sys.path.insert(0,'/home/great_czz/project/演化算法助教工作/算法评分赛/rscore')
from sector_math_check import Tree, score_brute, detach, attach

def parse_into_Tree(line):
    t=Tree(); lab={}; stack=[]; edges=[]
    i=0; n=len(line)
    while i<n:
        c=line[i]
        if c=='(':
            u=t.new(); stack.append(u); i+=1
        elif c in ',; \t\r\n': i+=1
        elif c==')':
            u=stack.pop(); i+=1; s=0
            if i<n and line[i] in 'SD':
                s=1 if line[i]=='S' else 0; i+=1
            lab[u]=s
            if stack: edges.append((stack[-1],u))
        elif c.isdigit():
            j=i
            while j<n and line[j].isdigit(): j+=1
            u=t.new(int(line[i:j])); i=j; lab[u]=1
            if stack: edges.append((stack[-1],u))
        else: i+=1
    for a,b in edges: t.edge(a,b,bool(lab.get(b,0)))
    return t

def py_truth(spfile, gfile):
    S=parse_into_Tree(open(spfile).read())
    genes=[parse_into_Tree(l) for l in open(gfile) if l.strip() and '(' in l]
    base=score_brute(S,genes)
    out={}
    for xn,xt in use_taxa(S):
        Sm,_=detach(S,xn)
        best=0
        for u in range(len(Sm.adj)):
            for v,_ in Sm.adj[u]:
                if u<v:
                    sc=score_brute(attach(Sm,xt,u,v),genes)
                    if sc-base>best: best=sc-base
        out[xt]=best
    return out

def use_taxa(S):
    return sorted(S.tax.items(), key=lambda kv: kv[1])

BIN="/home/great_czz/project/演化算法助教工作/算法评分赛/rscore/solver/sector_exact"
CASES=["c1_dominant","c2_all_small","c3_half","c4_star","c5_unary","c6_caterpillar","c7_missing","c8_nodash"]
tot=bad=0
for c in CASES:
    sp=f"/tmp/adv/{c}.sp.nwk"; g=f"/tmp/adv/{c}.g.trees"
    truth=py_truth(sp,g)
    r=subprocess.run([BIN,"-t",sp,"-g",g,"--all","--genes","100","-j","2","--scorer","/home/great_czz/project/演化算法助教工作/算法评分赛/rscore/solver/score_fast"],
                     capture_output=True,text=True)
    got={}
    for line in r.stdout.strip().split("\n"):
        f=line.split("\t")
        if len(f)==4: got[int(f[0])]=int(f[1])
    if not got:
        print(f"{c}: sector_exact FAILED: {r.stderr.strip()[:120]}"); bad+=1; continue
    diffs=[(t,truth[t],got.get(t)) for t in truth if got.get(t)!=truth.get(t)]
    tot+=len(truth); bad+=len(diffs)
    tag="OK " if not diffs else "DIFF"
    print(f"{c:16s} {tag} leaves={len(truth)}  " +
          ("" if not diffs else f"{len(diffs)} mismatch, e.g. {diffs[:3]}"))
print(f"=== {tot-bad}/{tot} (leaf, bestDelta) pairs agree ===")
