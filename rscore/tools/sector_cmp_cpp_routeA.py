import sys
sys.path.insert(0,'/home/great_czz/project/演化算法助教工作/算法评分赛/rscore')
from sector_math_check import (Tree, solid_cluster, all_dirs, solid_visible_from)

def T_formula_gamma(g, xnode, Ds, d, gamma):
    """Theorem 1 with A_x replaced by an explicit cluster `gamma` (route A)"""
    d2, d3 = [i for i in range(3) if i != d]
    total = 0
    for q in range(len(g.adj)):
        if q in g.tax: continue
        dirs = all_dirs(g, q)
        i0 = None
        for idx,(w,tx) in enumerate(dirs):
            if g.tax[xnode] in tx: i0 = idx; break
        if i0 is None: continue
        P = len(Ds[d] & dirs[i0][1] & gamma)
        if P == 0: continue
        q2,q3 = [],[]
        for idx,(w,tx) in enumerate(dirs):
            if idx == i0: continue
            vis = solid_visible_from(g,q,w)
            q2.append(len(Ds[d2] & tx & vis)); q3.append(len(Ds[d3] & tx & vis))
        cps = sum(q2[j]*q3[l] + q3[j]*q2[l] for j in range(len(q2)) for l in range(j+1,len(q2)))
        total += P*cps
    return total

def parse_tree_line(line):
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

def largest_cluster(g):
    best=None
    for u in g.tax:
        c=solid_cluster(g,u)
        if best is None or len(c)>len(best): best=c
    return best

S = parse_tree_line(open('/tmp/small_sp.nwk').read())
genes=[parse_tree_line(l) for l in open('/tmp/small_g.trees') if l.strip() and '(' in l]
taxa=sorted(S.tax.values())
ch=True
while ch:
    ch=False
    for u in range(len(S.adj)):
        if len(S.adj[u])==2:
            a,b=S.adj[u][0][0],S.adj[u][1][0]
            S.adj[a]=[(v,ss) for v,ss in S.adj[a] if v!=u]
            S.adj[b]=[(v,ss) for v,ss in S.adj[b] if v!=u]
            S.edge(a,b); S.adj[u]=[]; ch=True
xt=taxa[0]; xn=[u for u,t in S.tax.items() if t==xt][0]; p=S.adj[xn][0][0]
dump={}
for line in open('/tmp/small_dump.txt'):
    f=line.split()
    if f[0]=='T': dump[(int(f[1]),int(f[2]))]=(int(f[3]),int(f[4]),int(f[5]))
gam=[largest_cluster(g) for g in genes]
nbad=ncmp=0
for m in sorted(set(mm for mm,_ in dump)):
    if m==p: continue
    Ds=[d[1] for d in all_dirs(S,m)]
    for x in taxa:
        xn2=[u for u,t in S.tax.items() if t==x][0]
        DsT=[frozenset(s-{x}) for s in Ds]
        tri=[]
        for d in range(3):
            v=0
            for gi,g in enumerate(genes):
                gn=g.t2n()
                if x not in gn: continue
                v+=T_formula_gamma(g, gn[x], DsT, d, gam[gi])
            tri.append(v)
        got=dump[(m,x)]; ncmp+=1
        if sorted(tri)!=sorted(got):
            nbad+=1
            if nbad<=3: print(f"  MISMATCH m={m} x={x} py={sorted(tri)} cpp={sorted(got)}")
print(f"=== route-A(largest cluster) vs C++: {ncmp} pairs, {nbad} mismatches ===")
