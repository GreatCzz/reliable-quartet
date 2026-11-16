import sys, itertools
sys.path.insert(0,'/home/great_czz/project/演化算法助教工作/算法评分赛/rscore')
from sector_math_check import Tree, T_formula, solid_cluster, all_dirs

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

S = parse_tree_line(open('/tmp/small_sp.nwk').read())
genes=[parse_tree_line(l) for l in open('/tmp/small_g.trees') if l.strip() and '(' in l]
taxa = sorted(S.tax.values())
xt = taxa[0]
# suppress degree-2 nodes exactly like the C++ (keep ids)
ch=True
while ch:
    ch=False
    for u in range(len(S.adj)):
        if len(S.adj[u])==2:
            a,b=S.adj[u][0][0], S.adj[u][1][0]
            S.adj[a]=[(v,ss) for v,ss in S.adj[a] if v!=u]
            S.adj[b]=[(v,ss) for v,ss in S.adj[b] if v!=u]
            S.edge(a,b); S.adj[u]=[]; ch=True
xn=[u for u,t in S.tax.items() if t==xt][0]
p=S.adj[xn][0][0]
print("species leaves", len(S.tax), "internal", len([u for u in range(len(S.adj)) if u not in S.tax and S.adj[u]]))
print("x node", xn, "parent", p)

# read the C++ dump
dump={}
cur=None
for line in open('/tmp/small_dump.txt'):
    f=line.split()
    if f[0]=='M': cur=int(f[1])
    elif f[0]=='T':
        m=int(f[1]); t=int(f[2]); vals=(int(f[3]),int(f[4]),int(f[5]))
        dump[(m,t)]=vals
print("dump m's:", sorted(set(m for m,_ in dump)))

nbad=0; ncmp=0
for m in sorted(set(mm for mm,_ in dump)):
    if m==p: 
        print(f"  skip m={m} (== parent of x)"); continue
    dirs=all_dirs(S,m)
    Ds=[d[1] for d in dirs]
    for x in taxa:
        xn2=[u for u,t in S.tax.items() if t==x][0]
        # true tripartition of (S-x)-m : drop x from its component
        DsT=[frozenset(s-{x}) for s in Ds]
        tri=[]
        for d in range(3):
            v=0
            for g in genes:
                gn=g.t2n()
                if x not in gn: continue
                v+=T_formula(g, gn[x], DsT, d)
            tri.append(v)
        got=dump[(m,x)]
        ncmp+=1
        if sorted(tri)!=sorted(got):
            nbad+=1
            if nbad<=5:
                print(f"  MISMATCH m={m} x={x} py={sorted(tri)} cpp={sorted(got)}")
print(f"=== compared {ncmp} (m,x) pairs, {nbad} mismatches ===")
