"""Compare sector_sweep's T table against the exact per-leaf formula (Python)."""
import sys, subprocess
ADV="/home/great_czz/project/演化算法助教工作/算法评分赛"
sys.path.insert(0,ADV+"/rscore")
from sector_math_check import Tree, T_formula, all_dirs
SW=ADV+"/rscore/solver/sector_sweep"
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

def run(case):
    sp=f"{D}/{case}.sp.nwk"; g=f"{D}/{case}.g.trees"
    S=suppress_deg2(parse_into_Tree(open(sp).read()))
    genes=[parse_into_Tree(l) for l in open(g) if l.strip() and '(' in l]
    dump="/tmp/w/T.txt"
    r=subprocess.run([SW,"--species",sp,"--genes",g,"--strategy","m","--threads","1",
                      "--m-count","100000","--dump",dump],capture_output=True,text=True)
    if r.returncode!=0:
        print(case,"sector_sweep failed:",r.stderr.strip()[:200]); return None
    T={}
    for line in open(dump):
        f=line.split()
        if f and f[0]=='T':
            T[(int(f[1]),int(f[2]))]=(int(f[3]),int(f[4]),int(f[5]))
    ms=sorted({m for m,_ in T})
    taxa=sorted(S.tax.values())
    bad=0; tot=0; badleaves={}
    for m in ms:
        if len(S.adj[m])!=3: continue
        Ds=[d[1] for d in all_dirs(S,m)]
        for x in taxa:
            xn=[u for u,t in S.tax.items() if t==x][0]
            p=S.adj[xn][0][0]
            if m==p: continue            # S-x has no such node
            DsT=[frozenset(s-{x}) for s in Ds]
            tri=[]
            for d in range(3):
                v=0
                for gg in genes:
                    gn=gg.t2n()
                    if x not in gn: continue
                    v+=T_formula(gg, gn[x], DsT, d)
                tri.append(v)
            got=T.get((m,x))
            tot+=1
            if got is None or sorted(tri)!=sorted(got):
                bad+=1; badleaves[x]=badleaves.get(x,0)+1
    return tot,bad,badleaves

for case in ["c1_dominant","c2_all_small","c3_half","c4_star","c5_unary","c6_caterpillar","c7_missing","c8_nodash"]:
    r=run(case)
    if r:
        tot,bad,bl=r
        print(f"{case:16s} {tot-bad}/{tot} (m,x) agree  " + ("" if not bad else f"bad leaves: {dict(sorted(bl.items())[:6])}"))
