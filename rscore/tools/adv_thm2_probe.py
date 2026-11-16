"""Is Delta = SUM_m [T(m,dir(R)) - T(m,dir(R0))] true when T is taken from the DEFINITION?"""
import sys
ADV="/home/great_czz/project/演化算法助教工作/算法评分赛"
sys.path.insert(0,ADV+"/rscore")
from sector_math_check import Tree, T_direct, all_dirs
D="/tmp/w/med"
exec(open('/tmp/w/pin.py').read().split('S=suppress_deg2')[0].split('from sector_math_check')[1].split('\n',1)[1])
S=suppress_deg2(parse_into_Tree(open(f"{D}/m20_dom_1.sp.nwk").read()))
genes=[parse_into_Tree(l) for l in open(f"{D}/m20_dom_1.g.trees") if l.strip() and '(' in l]
gtas=[frozenset(g.tax.values()) for g in genes]
ms=[m for m in range(len(S.adj)) if m not in S.tax and len(S.adj[m])==3]
x=0
xn=[u for u,t in S.tax.items() if t==x][0]
p=S.nbrs(xn)[0]
def dirOf(m,uu,vv):
    for i,w in enumerate(S.nbrs(m)):
        if w==uu or w==vv: return i
    return -1
def Tsum(uu,vv):
    tot=0
    for m in ms:
        if m==p: continue
        d=dirOf(m,uu,vv)
        if d<0: continue
        dirs=all_dirs(S,m)
        Ds=[frozenset(t) for _,t in dirs]
        v=0
        for gi,g in enumerate(genes):
            if x not in g.t2n(): continue
            v+=T_direct(g, g.t2n()[x], Ds, d)
        tot+=v
    return tot
def Tsum_from_pretty(uu,vv):
    """same but the pair sets are (S-x) components (x removed)"""
    tot=0
    for m in ms:
        if m==p: continue
        d=dirOf(m,uu,vv)
        if d<0: continue
        dirs=all_dirs(S,m)
        Ds=[frozenset(t)-{x} for _,t in dirs]
        v=0
        for gi,g in enumerate(genes):
            if x not in g.t2n(): continue
            v+=T_direct(g, g.t2n()[x], Ds, d)
        tot+=v
    return tot
ab=[w for w in S.nbrs(p) if w!=xn]
print("leaf0 node",xn,"parent",p,"ab",ab)
for mode,fn in [("Ds with x",Tsum),("Ds without x",Tsum_from_pretty)]:
    a=fn(p,ab[0]); b=fn(12,14)
    print(f"  {mode}: T(R0)={a} T(R)={b} delta={b-a}   (exact delta = 32)")
