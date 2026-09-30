# 2026-W40 工作记录（09-28 至 10-04）

## 主题

完成二维非线性四分之一圆环的 membrane 与 bending 一致切线，并接入载荷步和
Newton–Raphson 迭代。

## 本周完成

- [x] 当前构型一阶、二阶导数的 RKPM 离散；
- [x] nonlinear membrane 内力、材料刚度和几何刚度；
- [x] 三维 binormal/normal 几何矩阵与二维位移自由度组合；
- [x] bending strain、内力和材料刚度；
- [x] 以局部自由度指标 `p` 构造 `delta_A`、`delta_B1/2`、`delta_C`、
  `delta_D1/2` 和 `delta_Bb`；
- [x] bending 几何刚度逐列组装；
- [x] 全局内力、切线刚度、约束和载荷步 Newton 迭代；
- [x] MATLAB 完整运行检查。

## 当前运行结果

- 20 个载荷步均收敛；
- 每个载荷步 4 次 Newton 迭代；
- 最终自由自由度残量：约 `4.14e-11`；
- 最终约束残量：约 `6.14e-18`；
- 总切线对称性误差：约 `7.06e-17`；
- 当前 nonlinear 载荷点位移：约 `-7.68375e-03`。

## 代码快照

- `code/matlab_linear_reference/main_quarter_ring_nonlinear_learning.m`
- `code/matlab_linear_reference/main_quarter_ring_learning.m`
- `code/matlab_linear_reference/rkpm_shape_1d.m`
- `code/matlab_linear_reference/NONLINEAR_MIGRATION_GUIDE.md`

## 下一步

- [ ] 使用中心差分验证完整 bending tangent；
- [ ] 检查小载荷极限与线性基准的一致性；
- [ ] 整理局部张量循环的理论笔记；
- [ ] 开始 Raknes 4.2 端弯矩大转动算例。
