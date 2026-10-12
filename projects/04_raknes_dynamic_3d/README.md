# 04 Raknes 三维动力学

## 目标算例

1. Raknes 4.3.2：三根铰接三维重力缆索。
2. Raknes 4.3.4：预张力双曲抛物面缆网。

## 模型范围

- 每节点三个位移自由度。
- 膜作用和密切平面内的标量弯曲。
- 不含独立扭转自由度。
- 静力切线可能奇异；检查时间离散后的有效矩阵。

## 当前已完成的动力学演示

`code/matlab_raknes_dynamic/main_helix_raknes_dynamic.m` 是从
`projects/05_full_spatial_rod` 中的自写 3D 螺旋线 rank-probe 代码复制出的动力学版本，
使用 Newmark 平均加速度格式（`beta=0.25`、`gamma=0.5`）。每个时间步内仍然使用
Newton 迭代，动力学有效切线为

\[
K_{\mathrm{eff}}
=K_t+\frac{M}{\beta\Delta t^2}
 +\frac{\gamma C_d}{\beta\Delta t}.
\]

当前演示采用小载荷 `totalForce = 1.0e-6`、40 个时间步和质量比例阻尼，
目的是先验证动力学质量项是否消除静力切线的秩亏。运行后会生成：

- `results/helix_dynamic.avi`：螺旋线变形动画；
- `results/helix_dynamic_tip_history.csv`：加载点位移随时间的历史；
- `results/helix_dynamic_result.mat`：位移、速度、加速度、质量矩阵、切线和秩诊断数据。

本次运行的诊断为：静力约束切线秩 `268/291`，动力学有效切线秩 `291/291`；
最终自由度平衡残差约为 `1.27e-12`，约束残差约为 `1.66e-22`。

## 必做诊断

- 质量矩阵与有效切线矩阵的秩。
- 时间步敏感性。
- 有/无弯曲稳定的对照。
- 动力松弛稳态与论文结果的对照。
