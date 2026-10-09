# RKPM 曲线结构研究路线

这七个目录表示理论和程序能力的递进关系，不是七套互不相关的程序。

```text
01_linear_2d
02_nonlinear_static_2d        <- 当前阶段
03_nonlinear_dynamic_2d
04_raknes_dynamic_3d
05_full_spatial_rod
06_multirod_contact
07_friction
```

## 总原则

每一阶段都按以下顺序推进：

```text
连续理论 -> 离散公式 -> 单点检查 -> 组装 -> 基准算例 -> 收敛/能量检查
```

只有当前阶段通过验收，才把它复制到 `weekly/年份/周次/` 作为稳定快照，然后进入下一阶段。

## 共用参考

- Raknes：位移型、无扭转、密切平面弯曲的非线性缆索。
- Wang–Bazilevs：PCA、RKPM、高阶导数与稳定节点积分。
- Greco–Cuomo：三个位移加独立扭转角的空间 Kirchhoff–Love 杆。
- Armero：director 形式的完整三维非线性 Kirchhoff 杆。
- Karapiperis：多梁线对线接触与摩擦。
