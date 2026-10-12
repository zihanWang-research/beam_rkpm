# 03 二维非线性动力学

## 目标算例

1. Raknes 4.3.1：二维重力简支缆索。
2. Raknes 4.3.3：二维铰接摆动缆索。

## 新增能力

- 一致或集总质量矩阵。
- 速度、加速度与初始条件。
- generalized-alpha/Newmark 时间积分。
- 质量比例阻尼和动力松弛。
- 位移、速度、能量的时间历程检查。

## 当前起点

你自己写的七百多行单文件 nonlinear 学习代码已经复制到：

`code/matlab_nonlinear_dynamic/`

该目录目前只有主程序和一个 RKPM 形函数文件：

- `main_quarter_ring_nonlinear_learning.m`（788 行）；
- `rkpm_shape_1d.m`。

它是静力版本的独立副本，尚未加入质量矩阵或时间积分，因此不会影响
`projects/02_nonlinear_static_2d/` 中已经验证的代码。后续将按照以下顺序逐步修改：

1. 保留现有 membrane、bending、内力和一致切线函数；
2. 增加二维一致质量矩阵；
3. 先实现线性 Newmark 时间步；
4. 再把非线性内力和 Newton 迭代放入每个时间步；
5. 最后对照 Raknes 4.3.1 的二维动力学算例。

## 进入条件

`02_nonlinear_static_2d` 的内力和一致切线已经通过验证。
