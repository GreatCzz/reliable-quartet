"""Ground truth (official ref_brute, exhaustive over positions) vs sector_exact."""
import sys, subprocess, re
ADV="/home/great_czz/project/演化算法助教工作/算法评分赛"
sys.path.insert(0,ADV+"/scorer")
sys.path.insert(0,ADV+"/rscore")
import ref_brute as RB
from sector_math_check import Tree, detach, attach
BIN=ADV+"/rscore/solver/sector_exact"
SF=ADV+"/rscore/solver/score_fast"
D=ADV+"/rscore/tools/adv_cases"

def parse_into_Tree(line):
    t=Tree(); lab={}; stack=[]; edges=[]
    i=0; n=len(line)
    while i<n:
        c=line[i]
        if c=='(': u=t.new(); stack.append(u); i+=1
        elif c in ',; \t\r\n': i+=1
        elif c==')':
            u=stack.pop(); i+=1; s=0
            if i<n and line[i] in 'SD': s=1 if line[i]=='S' else 0; i+=1
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

def nw(t):
    internals=[u for u in range(len(t.adj)) if u not in t.tax and len(t.adj[u])>=2]
    root=internals[0]
    out=[]
    def rec(u,p):
        if u in t.tax: out.append(str(t.tax[u])); return
        out.append('('); first=True
        for v,s in t.adj[u]:
            if v==p: continue
            if not first: out.append(',')
            first=False; rec(v,u)
            if v not in t.tax: out.append('S' if s else 'D')
        out.append(')')
    rec(root,-1); out.append(';'); return ''.join(out)

def genes_of(path):
    return [parse_into_Tree(l) for l in open(path) if l.strip() and '(' in l]

def score(S,genes):
    sp=nw(S)
    return sum(RB.brute_score(sp, nw(g))[0] for g in genes)

def suppress_deg2(t):
    ch=True
    while ch:
        ch=False
        for u in range(len(t.adj)):
            if len(t.adj[u])==2:
                a,b=t.adj[u][0][0],t.adj[u][1][0]
                t.adj[a]=[(v,ss) for v,ss in t.adj[a] if v!=u]
                t.adj[b]=[(v,ss) for v,ss in t.adj[b] if v!=u]
                t.edge(a,b); t.adj[u]=[]; ch=True
    return t

def truth(spfile,gfile):
    S=suppress_deg2(parse_into_Tree(open(spfile).read())); genes=genes_of(gfile)
    base=score(S,genes); out={}
    for xn,xt in sorted(S.tax.items(), key=lambda kv: kv[1]):
        Sm,_=detach(S,xn); best=0
        for u in range(len(Sm.adj)):
            for v,_ in Sm.adj[u]:
                if u<v:
                    d=score(attach(Sm,xt,u,v),genes)-base
                    if d>best: best=d
        out[xt]=best
    return out, base

CASES=["c1_dominant","c2_all_small","c3_half","c4_star","c5_unary","c6_caterpillar","c7_missing","c8_nodash"]
tot=bad=0
for c in CASES:
    sp=f"{D}/{c}.sp.nwk"; g=f"{D}/{c}.g.trees"
    tr,pybase=truth(sp,g)
    out=subprocess.run([BIN,"-t",sp,"-g",g,"--all","--genes","100","-j","2","--scorer",SF],
                       capture_output=True,text=True).stdout
    got={}
    for line in out.strip().split("\n"):
        f=line.split("\t")
        if len(f)==4: got[int(f[0])]=int(f[1])
    sfbase=None
    r=subprocess.run([SF,"-t",sp,"-g",g,"-j","1"],capture_output=True,text=True).stdout
    m=re.search(r'"score":(\d+)',r); sfbase=int(m.group(1)) if m else None
    diffs=[(t,tr[t],got.get(t)) for t in tr if got.get(t)!=tr[t]]
    tot+=len(tr); bad+=len(diffs)
    print(f"{c:16s} base ref={pybase} score_fast={sfbase} {'OK' if pybase==sfbase else 'BASE-MISMATCH'} | "
          f"leaves={len(tr)} {'OK' if not diffs else f'{len(diffs)} diff e.g. {diffs[:2]}'}")
print(f"=== {tot-bad}/{tot} (leaf,bestDelta) pairs agree ===")
