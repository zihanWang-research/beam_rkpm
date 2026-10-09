# 01 二维线性基准（已完成）

## 已验证内容

- Raknes 4.1.1、4.1.2、4.1.3 节点积分算例。
- 4.1.1、4.1.2 Gauss 积分。
- RKPM 阶次、影响域和 L2 收敛率扫描。
- 节点积分膜稳定项与 SCNI 试验。

## 主代码入口

- 节点积分：`../../code/matlab/benchmarks/raknes_4_1_linear/`
- Gauss 积分：`../../weekly/2026/2026-W38_09-14_to_09-20/code/matlab_gauss/`
- 稳定节点积分：`../../weekly/2026/2026-W38_09-14_to_09-20/code/matlab_stabilized_nodal/`
- 四分之一圆环个人学习版：`../../code/matlab/benchmarks/quarter_ring/user_1d_learning/`

## 回归测试职责

以后修改非线性、动力学或空间杆程序时，二维线性极限必须重新退化到本阶段结果。
