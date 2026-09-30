# 二维非线性静力周快照

本目录记录 2026-W39 开始的 nonlinear 学习算例。

对应的小变形起点副本位于同级目录 `../matlab_linear_reference/`；后续学习修改
从该副本开始，原始线性 benchmark 保持不变。

- `main_nonlinear_membrane_bar_learning.m`：直杆轴向拉伸初始化骨架。
- `main_nonlinear_quarter_ring_learning.m`：已同步的完整四分之一圆环 nonlinear 基线，包含 membrane、planar bending、分步加载和 Newton 迭代。
- 活跃代码位于 `../../../../../code/matlab/nonlinear/static_2d/`。

本周快照由活跃代码同步生成；后续实现优先写入活跃代码，通过阶段检查后再更新周记录。

学习顺序固定为：直杆 membrane → 四分之一圆环 membrane/bending → Raknes 4.2。
