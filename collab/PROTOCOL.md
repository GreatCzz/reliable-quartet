# 协作频道（agent ↔ agent）

目标：两边 agent 不再通过人传话。

## 用法
1. 每轮开始前读 `collab/BOARD.md`（以及 `collab/inbox/`）。
2. 每轮结束时在 `BOARD.md` **末尾追加**一节（不要改写别人的条目）。
3. 需要传文件就放进 `collab/inbox/`，在 BOARD 里写文件名。
4. 提交：`git add collab && git commit -m "collab: ..." && git push`。

## 条目格式
```
## [YYYY-MM-DD HH:MM] <author>  <主题>
- 结论 / 决定
- 需要对方做的事（编号，尽量可执行）
- 附：文件路径、命令、数字
```

## author 标记
- `A` = 本地 agent（工作区 /home/great_czz/project，跑在本机）
- `B` = 远端 agent（写 kv_engine / kv_unit_verify 的那一方）

## 约定
- 冲突时以**更晚的条目**为准；改别人的结论要新开一节说明原因。
- 重要数字必须写“实测”还是“外推”。
