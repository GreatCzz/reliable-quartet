# 环境建议

## 模型
- 当前会话默认 `deepseek-flash`；本项目可用 **`deepseek-v4-pro`**（更强，贵）；
- 建议：**算法设计/核心 C++/对拍**用 v4-pro；跑数据、写文档用 flash；
- subagent 也可以单独指定 v4-pro（并行做独立任务）。

## 会话
- 每个子任务**新开会话**，先读本文件夹 + `AGENTS.md`；
- 不要把跑数据、等待、调试混在一个超长会话里（>50 轮质量明显下降）。

## 计算
- 本地：32 核；全量 NNI 一轮 1 小时，建议后台跑；
- 集群：`aster.pku.edu.cn`（`doc/cluster/aster.md`），SLURM，适合并行多参数/多起点；
- 跑大任务前确认电源（电池模式下会降频，实测 3 倍差距）。

## 工具
- 固定基准：`handoff/bench.sh`（500 类群 + 1000 类群）；
- 回归：`./verify.sh`；
- 性能剖析：`perf record -g ./solver ...`。
