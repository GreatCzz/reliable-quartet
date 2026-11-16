# 复现：官方暴力参考 vs score_fast
cd "/home/great_czz/project/演化算法助教工作/算法评分赛"
python3 - <<'PY'
import sys, subprocess, re
sys.path.insert(0,'scorer'); import ref_brute as R
for c in ["c1_dominant","c3_half","c4_star","c5_unary","c6_caterpillar","c8_nodash"]:
    sp=open(f"rscore/tools/adv_cases/{c}.sp.nwk").read().strip()
    gs=[l.strip() for l in open(f"rscore/tools/adv_cases/{c}.g.trees") if l.strip()]
    ref=sum(R.brute_score(sp,g)[0] for g in gs)
    out=subprocess.run(["rscore/solver/score_fast","-t",f"rscore/tools/adv_cases/{c}.sp.nwk","-g",f"rscore/tools/adv_cases/{c}.g.trees","-j","1"],capture_output=True,text=True).stdout
    print(f"{c:16s} ref_brute={ref:<6d} score_fast={re.search(chr(34)+'score'+chr(34)+r':(\\d+)',out).group(1)}")
PY
