# 代码架构与构建

## 目录（上一级 = 项目根 `算法评分赛/`）

```
scorer/                  Python 参考实现（暴力 + 快速公式 + 对拍）—— 唯一权威口径
rscore/
  solver/
    solver.cpp           ★ 我们的 NNI 局部搜索（W 正确）+ 打分
    score.cpp            参考打分器（慢，6 分钟/全量）
    score_fast.cpp       快速打分器（导师算法移植，5 秒/全量）★ 默认
    builder.cpp          打分插入构造（有 bug 历史，已被队友修过一版）
    recombine.cpp        三分划并集 DP 重组
    incr_dp.cpp          增量 13 状态 DP 原型（已验证）
  aster2/                ASTER2 工作副本 + rscore.hpp（O(n) class）+ rscore_incr.hpp（增量）
  aster1/                ASTER v1 源码（只读，ASTRAL-IV 可编译）
  score.sh               只评分入口（默认 fast，--ref / --verify）
  rbuild.sh              一键：共识 → NNI → 打分
verify.sh                全仓自检
submissions/             提交树 + SCOREBOARD.md
启发式算法评分赛全体数据/   rep01..rep20 基因树（gz）
```

## 构建

```bash
cd rscore/solver
g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src solver.cpp     -o solver
g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src score.cpp      -o score
g++ -std=c++20 -O3 -march=native -fopenmp score_fast.cpp                  -o score_fast
g++ -std=c++20 -O3 -march=native -fopenmp incr_dp.cpp                     -o incr_dp
cd ../aster2 && make rscore_incr
cd ../../aster1 && make astral        # ASTRAL-IV（需先去掉单叉节点）
```

## 运行

```bash
# 只评分（默认快后端，JSON）
./rscore/score.sh -t tree.nwk -g genes.trees.gz -j 32

# 我们的 NNI
./rscore/solver/solver -s tree.nwk -i genes.trees -t 32 --sweeps 2 -o out.nwk

# ASTER2 框架 + 增量计分（仅放置）
./rscore/aster2/bin/rscore_incr -i genes.trees -t 16 -o out.nwk

# ASTRAL-IV（先 splice 单叉节点）
cd aster1 && ./bin/astral4 -i genes_spliced.trees -o out.tree -t 32
```
