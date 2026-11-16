import random, sys
sys.path.insert(0,'/home/great_czz/project/演化算法助教工作/算法评分赛/rscore')
from sector_math_check import rand_binary, rand_gene

def to_newick_binary(S):
    # unrooted binary array -> newick string rooted at the first internal node
    # S: Tree with adj and tax
    internal = [u for u in range(len(S.adj)) if u not in S.tax]
    root = internal[0]
    out=[]
    def rec(u, p):
        if u in S.tax:
            out.append(str(S.tax[u])); return
        out.append('(')
        first=True
        for v,_ in S.adj[u]:
            if v==p: continue
            if not first: out.append(',')
            first=False
            rec(v,u)
        out.append(')')
    rec(root,-1)
    out.append(';')
    return ''.join(out)

def gene_newick(g):
    # g is an unrooted Tree (Tree with adj+tax) -> newick rooted at internal[0]
    internal=[u for u in range(len(g.adj)) if u not in g.tax]
    root=internal[0]
    out=[]
    def rec(u,p):
        if u in g.tax:
            out.append(str(g.tax[u])); return
        out.append('(')
        first=True
        for v,s in g.adj[u]:
            if v==p: continue
            if not first: out.append(',')
            first=False
            rec(v,u)
            if s: out.append('S')   # label goes right after the child's ')'
        out.append(')')
    rec(root,-1)
    out.append(';')
    return ''.join(out)

rng=random.Random(int(sys.argv[1]))
n=9
S=rand_binary(n,rng)
genes=[]
for _ in range(3):
    taxa=[a for a in range(n) if rng.random()>0.1]
    if len(taxa)<5: taxa=list(range(5))
    genes.append(rand_gene(taxa,rng,float(sys.argv[2])))
open('/tmp/small_sp.nwk','w').write(to_newick_binary(S)+'\n')
with open('/tmp/small_g.trees','w') as f:
    for g in genes: f.write(gene_newick(g)+'\n')
print("wrote")
