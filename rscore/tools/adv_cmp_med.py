"""Medium-scale: sector_exact (truth) vs sector_place (fast) best-delta per leaf."""
import subprocess, glob, os, re
ADV="/home/great_czz/project/演化算法助教工作/算法评分赛"
EX=ADV+"/rscore/solver/sector_exact"; PL=ADV+"/rscore/solver/sector_place"
SW=ADV+"/rscore/solver/sector_sweep"; SF=ADV+"/rscore/solver/score_fast"
for sp in sorted(glob.glob("/tmp/w/med/*.sp.nwk")):
    name=os.path.basename(sp)[:-7]
    g=f"/tmp/w/med/{name}.g.trees"
    # truth
    r=subprocess.run([EX,"-t",sp,"-g",g,"--all","--genes","100","-j","4","--scorer",SF],
                     capture_output=True,text=True)
    truth={}
    for line in r.stdout.strip().split("\n"):
        f=line.split("\t")
        if len(f)==4: truth[int(f[0])]=int(f[1])
    # fast
    d="/tmp/w/T2.txt"
    r2=subprocess.run([SW,"--species",sp,"--genes",g,"--strategy","m","--threads","2",
                       "--m-count","100000","--dump",d],capture_output=True,text=True)
    r3=subprocess.run([PL,"--species",sp,"--T",d,"--all","--apply","0"],capture_output=True,text=True)
    fast={}
    for line in r3.stdout.strip().split("\n"):
        f=line.split()
        if len(f)==3 and f[0]=="A": fast[int(f[1])]=int(f[2])
    if not truth:
        print(f"{name}: sector_exact produced nothing: {r.stderr.strip()[:120]}"); continue
    diffs=[(t,truth[t],fast.get(t)) for t in truth if fast.get(t)!=truth[t]]
    print(f"{name:16s} leaves={len(truth)} fast_reported={len(fast)} mismatches={len(diffs)}"
          + ("" if not diffs else f"  e.g. {diffs[:3]}"))
