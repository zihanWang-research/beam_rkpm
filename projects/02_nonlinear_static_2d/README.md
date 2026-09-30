# 02 二维非线性静力（当前阶段）

## 当前状态（2026-09-30）

二维非线性四分之一圆环已经形成可运行版本：

- Green–Lagrange membrane strain、膜内力及材料/几何切线；
- 基于当前 binormal/normal 标架的 bending strain 和 bending 内力；
- bending 材料切线与一致几何切线；
- 载荷步、约束增广系统和 Newton–Raphson 迭代；
- 三维几何矩阵与二维位移自由度的组合实现。

当前主程序在 20 个载荷步中均于 4 次 Newton 迭代内收敛。最终自由自由度残量约为
`4.14e-11`，约束残量约为 `6.14e-18`，切线对称性误差约为 `7.06e-17`。

## 下一目标

1. 对 bending 一致切线执行有限差分验证；
2. 检查小载荷极限是否回到线性四分之一圆环结果；
3. 复现 Raknes 4.2：端部恒定弯矩作用下的二维悬臂梁大转动。

## 学习与实现顺序

1. 参考构型、当前构型和位移。
2. 非线性膜应变与曲率变化。
3. 应变的一阶变分。
4. 内力残量。
5. 内力残量的一致线性化和切线刚度。
6. 载荷步与 Newton–Raphson 迭代。
7. 刚体转动、有限差分切线和论文圆形解验证。

## 起始代码

从 `../01_linear_2d/` 指向的二维 RKPM 程序中选择一份已验证主程序；保留 RKPM 形函数、积分和约束模块，替换线性应变、刚度与一次求解部分。

## 当前代码入口

- 线性学习基准副本：`code/matlab_linear_reference/main_quarter_ring_learning.m`。
- 当前 nonlinear 主程序：`code/matlab_linear_reference/main_quarter_ring_nonlinear_learning.m`。
- 逐阶段修改指南：`code/matlab_linear_reference/NONLINEAR_MIGRATION_GUIDE.md`。
- 第一步：`../../code/matlab/nonlinear/static_2d/main_nonlinear_membrane_bar_learning.m`
- 第一个验证：二维直杆轴向拉伸，只验证 nonlinear membrane。
- 第二步：`../../code/matlab/nonlinear/static_2d/main_nonlinear_quarter_ring_learning.m`，学习初始曲率下的 nonlinear membrane 与 bending。
- 当前完整四分之一圆环副本：`code/matlab_nonlinear/main_nonlinear_quarter_ring_learning.m`。
- 第三步：在内力和一致切线通过有限差分后，实现端弯矩并复现 Raknes 4.2。

## 完成标准

- 残量随 Newton 迭代稳定下降。
- 一致切线通过有限差分验证。
- 刚体运动不产生内力。
- 端弯矩算例能够形成论文中的圆弧构型。

## 运行

在 MATLAB 中进入 `code/matlab_linear_reference/` 后运行：

```matlab
main_quarter_ring_nonlinear_learning
```

主程序依赖同目录下的 `rkpm_shape_1d.m`。
