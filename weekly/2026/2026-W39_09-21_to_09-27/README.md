# 2026-W39 工作计划（09-21 至 09-27）

## 主题：二维几何非线性与 nonlinear membrane

本周开始 `02_nonlinear_static_2d`。先用直杆完成膜内力和 Newton 闭环，再进入四分之一圆环和 Raknes 4.2。

## 必须完成

- [ ] 推导并实现 Green–Lagrange 膜应变、当前 `Bm` 和膜内力。
- [ ] 分别组装材料切线 `Kmat` 和几何切线 `Kgeo`。
- [ ] 用中心差分验证膜内力 Jacobian。
- [ ] 完成直杆载荷步和 Newton–Raphson。
- [ ] 与均匀伸长解 `P = EA*lambda*(lambda^2-1)/2` 对照。

## 有余力时完成

- [ ] 将已验证膜链放入四分之一圆环。
- [ ] 推导二维非线性曲率和 bending 残量。
- [ ] 对弯曲切线做有限差分检查。

## 下一阶段

Raknes 4.2 端部恒定弯矩大转动悬臂梁不是本周必须完成项，应在膜和曲率切线通过后再开始。

## 代码位置

- 活跃实现：`../../../code/matlab/nonlinear/static_2d/`
- 本周快照：`code/matlab_nonlinear/`

