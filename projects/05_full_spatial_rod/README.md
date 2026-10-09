# 05 完整空间杆

本阶段不直接覆盖 Raknes 程序，而是建立独立模型。

当前新增的秩测试分支位于：
`code/matlab_own_nonlinear_rank_probe/`。
该分支只沿用自己的二维 nonlinear membrane+bending Newton 框架，测试三维 PCA
螺旋线当前切线的秩，不使用导师的 Rodrigues/director 实现。

## 目标应变

```text
轴向应变 + 两个弯曲曲率 + 一个扭转应变
```

## 候选运动学

- Greco–Cuomo：三个位移加一个独立扭转角，较接近当前位移型 RKPM。
- Armero：中心线加截面 directors，更一般，但自由度和一致线性化更复杂。

## 验证要求

- 两个主方向纯弯曲。
- 直杆纯扭转。
- 叠加刚体转动的客观性。
- 三维弯曲与扭转耦合基准。

## 当前预检：3D 螺旋线切线秩

在进入完整 3D nonlinear rod 之前，已运行
`code/matlab/geometry/pca_nodal_integration/compare_helix_3d_bending_rank.m`。
当前离散下，约束后的自由子空间大小为 291：

| 模型 | 数值秩 | 零空间维数 |
|---|---:|---:|
| 膜内力 | 98 | 193 |
| 膜 + 标量 bending | 197 | 94 |
| 膜 + 两个法向 bending 方向 | 291 | 0 |

这说明原先只有一个标量弯曲应变时，3D 螺旋线仍然秩亏；加入法向平面内的两个弯曲方向后，当前诊断矩阵满秩。该结果是小变形秩预检，不等同于完整 3D nonlinear tangent 已经实现。
