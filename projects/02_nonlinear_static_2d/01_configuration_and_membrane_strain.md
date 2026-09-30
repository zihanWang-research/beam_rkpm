# 第一课：构型与非线性膜应变

## 本课位置

```text
连续理论
  -> 参考构型与当前构型          <- 1.1
  -> 非线性膜应变                <- 1.2
  -> 膜应变的一阶变分            <- 1.3
  -> 膜内力与一致切线            <- 1.4
  -> 弯曲应变及其线性化
  -> Newton-Raphson
```

建议学习时间：60～90 分钟。本课不写完整程序，只完成公式推导和几个标量检查。

## 1. 参考构型与当前构型

以材料参数 `xi` 标记中心线上的材料点。参考构型为

```math
\mathbf X=\mathbf X(\xi),
```

当前构型为

```math
\mathbf x=\mathbf x(\xi)=\mathbf X(\xi)+\mathbf u(\xi).
```

对参数求导：

```math
\mathbf x_{,\xi}
=
\mathbf X_{,\xi}+\mathbf u_{,\xi}.
```

参考和当前线元分别为

```math
\mathrm dS=\|\mathbf X_{,\xi}\|\,\mathrm d\xi,
\qquad
\mathrm ds=\|\mathbf x_{,\xi}\|\,\mathrm d\xi.
```

因此局部伸长比为

```math
\lambda
=
\frac{\mathrm ds}{\mathrm dS}
=
\frac{\|\mathbf x_{,\xi}\|}
     {\|\mathbf X_{,\xi}\|}.
```

## 2. Raknes 的非线性膜应变

Raknes 使用中心线 Green-Lagrange 应变的协变分量：

```math
\varepsilon
=
\frac12
\left(
\mathbf x_{,\xi}\cdot\mathbf x_{,\xi}
-
\mathbf X_{,\xi}\cdot\mathbf X_{,\xi}
\right).
```

代入 `x = X + u`：

```math
\begin{aligned}
\varepsilon
&=
\frac12
\left[
(\mathbf X_{,\xi}+\mathbf u_{,\xi})
\cdot
(\mathbf X_{,\xi}+\mathbf u_{,\xi})
-
\mathbf X_{,\xi}\cdot\mathbf X_{,\xi}
\right]\\
&=
\mathbf X_{,\xi}\cdot\mathbf u_{,\xi}
+
\frac12
\mathbf u_{,\xi}\cdot\mathbf u_{,\xi}.
\end{aligned}
```

其中

```math
\varepsilon_L
=
\mathbf X_{,\xi}\cdot\mathbf u_{,\xi}
```

是线性部分，而

```math
\varepsilon_{NL}
=
\frac12
\mathbf u_{,\xi}\cdot\mathbf u_{,\xi}
```

是几何非线性部分。以前的小变形程序忽略的正是第二项。

若参考参数就是参考弧长 `S`，则

```math
\|\mathbf X_{,S}\|=1,
\qquad
\varepsilon
=
\frac12(\lambda^2-1).
```

## 3. 膜应变的一阶变分

参考构型固定，因此

```math
\delta\mathbf X=\mathbf 0,
\qquad
\delta\mathbf x=\delta\mathbf u.
```

对膜应变取变分：

```math
\boxed{
\delta\varepsilon
=
\mathbf x_{,\xi}\cdot\delta\mathbf u_{,\xi}
}
```

与线性理论不同，这里的系数是当前构型导数 `x_,xi`，所以每次非线性迭代都要更新。

## 4. 膜内力与一致切线

若暂时以参考弧长为参数，膜内能为

```math
U_m
=
\frac12
\int_0^L EA\,\varepsilon^2\,\mathrm dS.
```

内虚功为

```math
\delta U_m
=
\int_0^L
EA\,\varepsilon\,\delta\varepsilon
\,\mathrm dS.
```

为了获得 Newton 法需要的切线，对 `delta U_m` 沿增量 `Delta u` 再线性化：

```math
\Delta\varepsilon
=
\mathbf x_{,S}\cdot\Delta\mathbf u_{,S},
```

```math
\Delta(\delta\varepsilon)
=
\delta\mathbf u_{,S}\cdot\Delta\mathbf u_{,S}.
```

所以

```math
\boxed{
\Delta(\delta U_m)
=
\int_0^L EA
\left[
(\delta\varepsilon)(\Delta\varepsilon)
+
\varepsilon
\,\delta\mathbf u_{,S}\cdot\Delta\mathbf u_{,S}
\right]
\mathrm dS
}
```

第一项产生材料刚度，第二项产生与当前轴向应力有关的几何刚度。

## 本课检查题

1. 当 `u = 0` 时，验证 `epsilon = 0`。
2. 对任意常量平移 `u = c`，验证 `epsilon = 0`。
3. 对刚体转动后的中心线，使用完整非线性应变验证 `epsilon = 0`。
4. 说明为什么只保留线性项时，大转动可能产生虚假应变。

完成这些检查以后，再进入弯曲应变 `j`、法向量变化 `delta n` 和弯曲切线刚度。
