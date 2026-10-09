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

局部 bending 切线和完整全局切线均已通过中心差分检查。完整全局切线的相对误差为
`1.57e-10`。解析式对应侧点 `A`（`theta=0`）的横向位移；总载荷为 `2` 时，
数值值为 `-7.00059530e-03`，解析值为 `-6.82848862e-03`，误差为 `2.52%`。
将总载荷降至 `2e-3` 后，同一点的误差为 `2.11%`，说明小载荷下 nonlinear
结果趋近线性解析解。之前使用顶部 `DY(end)` 与该解析式比较，属于点位不一致，
已经修正为侧点 `-DX(1)`。

在相同小载荷和相同比较点下，节点收敛结果为：

完整数据表：[nonlinear_node_convergence.csv](results/nonlinear_node_convergence.csv)

| 节点数 | 相对误差 |
|---:|---:|
| 17 | 2.11% |
| 33 | 1.1603% |
| 65 | 0.7159% |
| 129 | 0.5528% |

误差随节点加密持续下降，说明当前 nonlinear 实现已经通过基本的空间离散收敛检查；
剩余误差主要需要从影响域、积分和端点离散处理方面继续分析。

需要注意：`code/matlab/benchmarks/quarter_ring/user_1d_learning/node_integration_results/`
中的旧线性 CSV 是历史参数扫描结果，其中解析位移量级与当前 nonlinear 参数不一致，
不能直接用作线性极限对照。最终线性极限检查必须使用相同的
`E`、截面、`R/t`、总载荷、积分规则、边界条件和比较点重新运行。

同参数线性极限检查已经完成。线性程序在 `F=1` 时的位移斜率为
`3.38864014e-03`；nonlinear 程序在 `F=2e-6` 时的位移斜率为
`3.41960349e-03`，对应解析斜率 `3.41424431e-03`，相对误差为 `0.1570%`。
因此 nonlinear 在极小载荷下正确退化到线性响应。完整数据见
`results/linear_limit_check.csv`。

## 下一目标

1. 复现 Raknes 4.2：端部恒定弯矩作用下的二维悬臂梁大转动；
2. 整理最终代码和理论笔记。

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
- 当前 nonlinear 文件夹中的逐步完成版本：`code/matlab_nonlinear/main_quarter_ring_nonlinear_learning.m`。
- 局部 bending 切线检查：`code/matlab_nonlinear/checkQuarterRingBendingTangent.m`。
- 逐阶段修改指南：`code/matlab_linear_reference/NONLINEAR_MIGRATION_GUIDE.md`。
- 第一步：`../../code/matlab/nonlinear/static_2d/main_nonlinear_membrane_bar_learning.m`
- 第一个验证：二维直杆轴向拉伸，只验证 nonlinear membrane。
- 第二步：`../../code/matlab/nonlinear/static_2d/main_nonlinear_quarter_ring_learning.m`，学习初始曲率下的 nonlinear membrane 与 bending。
- 当前完整四分之一圆环副本：`code/matlab_nonlinear/main_nonlinear_quarter_ring_learning.m`。
- 第三步：在内力和一致切线通过有限差分后，实现端弯矩并复现 Raknes 4.2。

## 完成标准

- 残量随 Newton 迭代稳定下降。
- 局部 bending 和完整全局一致切线均通过有限差分验证。
- 刚体运动不产生内力。
- 端弯矩算例能够形成论文中的圆弧构型。

## 运行

在 MATLAB 中进入 `code/matlab_linear_reference/` 后运行：

```matlab
main_quarter_ring_nonlinear_learning
```

主程序依赖同目录下的 `rkpm_shape_1d.m`。
