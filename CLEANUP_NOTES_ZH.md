# 从原始工程中移除的内容

本清洁包的目标不是保存完整开发历史，而是提供可审计的论文复现材料。

## 删除类别

- Round-1/Round-2/V22/V23/V24A–V24J 等历史诊断与补丁脚本；
- Peng 独立复现的 Phase-2/Phase-3 开发文件；
- MESA qualification 的开发历史文件；
- 未进入最终主曲线的 third-party atomic-MMV 源码包；
- 旧版/重复绘图脚本，例如 `plot_fig4_resolution.m`、
  `plot_submission_final_core.m`、`plot_submission_final_core_v2.m`；
- 所有 diagnostic PNG、旧 MAT/CSV、quick/full 开发报告；
- Monte Carlo checkpoint 文件；
- 多轮 PATCH_SCOPE / SHA manifest / hotfix / install README；
- 离线 Python 审计脚本和仅用于开发决策的 CSV；
- 顶层历史运行器 `run_round*`, `run_v23*`, `run_v24*` 等；
- 早期 theoretical/Prony 辅助验证文件（不属于最终数值复现入口）。

## 为什么可以删

这些文件用于**方法开发、失败路径审计或阶段性资格冻结**，
但最终论文图和正式 paired Monte Carlo 不再调用它们。
删除它们不会改变清洁包中保存的最终数值结果，也不会改变最终绘图输入。

## 保留原则

只要满足下面任一条件就保留：
1. 最终主 Monte Carlo 直接调用；
2. 最终 8 张图直接调用或读取；
3. Fig.1/5/6/S1 的固定支持结果可由该脚本生成；
4. 支持实验运行前的资格门槛直接需要；
5. 最终论文报告的结果文件或 baseline qualification 记录。


## GitHub/Zenodo v2 additional reduction

Compared with the first clean package, this public-repository version also
removes the entire `plotting/` directory and the generated `figures/` folder.
The numerical MAT/CSV data underlying the figures remain included.

Reason: the public repository is focused on algorithm/data reproducibility,
while journal-specific publication layout is not required to verify the
reported numerical findings.
