# testdata

- `real40.tree` / `real40.trees`：rep01 裁到 40 类群（seed 7）+ 20 棵真实基因树（多叉、含缺失类群）

  ```bash
  /tmp/spr44_f rscore/testdata/real40.tree rscore/testdata/real40.trees 20   # 逐位置对拍，残差应为 0
  ```

> 早先的 `tree20d/g20d` 已删除：那批基因树是星形（无内部边），任何四叶树都不被显示，分数恒为 0，
> 逐位置比较退化成 `0 == 0`，是一份**假夹具**。
