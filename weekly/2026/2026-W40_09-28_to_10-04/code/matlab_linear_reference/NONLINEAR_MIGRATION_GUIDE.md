# 从线性四分之一圆环到二维几何非线性：修改路线

## 目标与文件规则

学习起点是本目录中的 `main_quarter_ring_learning.m`。这个文件作为线性基准，
后续不直接修改。正式学习时复制为：

`main_quarter_ring_nonlinear_learning.m`

完整的 `../main_nonlinear_quarter_ring_learning.m` 只作为最终参考，不在每一步开始前
照抄。每完成一个阶段，先通过本阶段检查，再进入下一阶段。

## 总体结构

线性程序的主线是：

`geometry -> RKPM -> Kmem/Kbend -> K*d=F -> plot`

非线性程序的主线要改成：

`geometry -> RKPM -> Fint(d) -> residual -> Kt(d) -> Newton -> plot`

真正改变的是力学核心和求解器。几何、RKPM、积分、自由度编号、约束矩阵和绘图框架
尽量保留。

---

## 阶段 0：冻结并验证线性基准

### 保留

- `nNode`、材料和截面参数；
- `thetaNode`、`X`、`h`；
- PCA/RKPM 邻域和形函数；
- `N`、`Neta`、`Netaeta`；
- `Kmem`、`Kbend`、`C`、`F`；
- 原来的位移恢复和绘图。

### 需要记录

- 线性程序是否正常运行；
- `uA_RKPM`；
- `uA_exact`；
- `relativeError_uA`；
- `norm(K-K','fro')/norm(K,'fro')`；
- `norm(C*d-g)`。

### 通过标准

程序无报错，位移、误差和矩阵对称性有可重复的基准值。

---

## 阶段 1：只整理程序结构，不改变线性结果

### 修改位置

原程序中组装 `Kmem`、`Kbend` 和求解 `Kaug\raug` 的部分。

### 修改内容

1. 把外力变量 `F` 改名为 `Fext`。
2. 保留原线性刚度：

   `Klinear = Kmem + Kbend`

3. 把线性内力明确写成：

   `Fint = Klinear*d`

4. 把平衡残量明确写成：

   `residual = Fint - Fext`

本阶段仍然只做一次线性求解，不加入 Newton。

### 为什么先做这一步

非线性程序的核心变量是 `Fint` 和 `residual`。先在线性问题中建立相同的软件结构，
可以把“程序重构错误”和“非线性公式错误”分开。

### 通过标准

重构前后的位移差满足：

`norm(dNew-dOld)/max(norm(dOld),eps) < 1e-12`

---

## 阶段 2：建立当前构型，但暂时不改 bending

### 新增运动学量

对每个积分点或节点积分点，建立局部自由度 `dLocal`，然后构造：

`uS  = H1*dLocal`

`uSS = H2*dLocal`

参考切向量和二阶导数记为：

`ARef = X_,S`

`CRef = X_,SS`

当前量为：

`aCurrent = ARef + uS`

`cCurrent = CRef + uSS`

对应公式：

\[
\mathbf a=\mathbf A+\mathbf u_{,S},\qquad
\mathbf c=\mathbf C+\mathbf u_{,SS}.
\]

### 本阶段只检查运动学

- `d=0` 时，`aCurrent=ARef`；
- `d=0` 时，`cCurrent=CRef`；
- 刚体平移时，`uS=0`、`uSS=0`；
- 参考切向量长度在弧长参数下应接近 1。

本阶段不组装新的内力。

---

## 阶段 3：把线性 membrane 替换为 nonlinear membrane

### 删除或停用

原来直接组装 `Kmem` 的代码块暂时保留在注释中用于对照，但不再用于总内力。

### 新增公式

膜应变：

\[
\mathcal E=\frac12(\mathbf a\cdot\mathbf a-1).
\]

膜应变的一阶变分：

\[
\delta\mathcal E=\mathbf a\cdot\delta\mathbf u_{,S}.
\]

离散形式：

\[
\delta\mathcal E=\mathbf B_m\delta\mathbf d,
\qquad
\mathbf B_m=\mathbf a^T\mathbf H_1.
\]

膜内力：

\[
\mathbf F_m^{\mathrm{int}}
=\int \mathbf B_m^T EA\mathcal E\,dS.
\]

### 过渡阶段的 bending

为了只检查 membrane，暂时保留线性弯曲内力：

`Fbend = Kbend*d`

总内力暂时写成：

`Fint = FmemNonlinear + Fbend`

### 通过标准

- `d=0` 时 `norm(FmemNonlinear)` 接近 0；
- 刚体平移不产生膜内力；
- 极小位移下 nonlinear membrane 内力与线性 membrane 内力接近；
- 局部和全局自由度编号没有维度错误。

---

## 阶段 4：加入 residual、数值切线和 Newton 迭代

### 残量

\[
\mathbf R(\mathbf d)
=\mathbf F^{\mathrm{int}}(\mathbf d)
-\lambda\mathbf F^{\mathrm{ext}}.
\]

代码结构：

`residual = Fint - loadFactor*Fext`

### 第一版切线

先使用中心差分：

\[
\mathbf K_t(:,i)\approx
\frac{\mathbf F^{\mathrm{int}}(\mathbf d+h\mathbf e_i)
-\mathbf F^{\mathrm{int}}(\mathbf d-h\mathbf e_i)}{2h}.
\]

### Newton 方程

\[
\begin{bmatrix}
\mathbf K_t&\mathbf C^T\\
\mathbf C&\mathbf0
\end{bmatrix}
\begin{bmatrix}
\Delta\mathbf d\\
\Delta\boldsymbol\lambda
\end{bmatrix}
=
\begin{bmatrix}
-\mathbf R\\
\mathbf g-\mathbf C\mathbf d
\end{bmatrix}.
\]

更新：

`d = d + deltaD`

### 通过标准

- 小载荷下 Newton 收敛；
- 每个载荷步的残量总体下降；
- `norm(C*d-g)` 足够小；
- 小载荷结果接近原线性结果。

---

## 阶段 5：加入载荷步

增加：

`for loadStep = 1:nLoadStep`

`loadFactor = loadStep/nLoadStep`

每个新载荷步使用上一个载荷步收敛后的 `d` 作为初值。

### 为什么需要

大变形问题的平衡路径是非线性的。一次施加完整载荷可能使 Newton 初值离解太远。

### 通过标准

- 所有载荷步都收敛；
- 减小载荷步长后最终解基本不变；
- 不依赖异常大的有限差分步长才能收敛。

---

## 阶段 6：把 linear bending 替换为 nonlinear bending

### 当前曲率

平面中心线曲率为：

\[
\kappa=
\frac{\det(\mathbf a,\mathbf c)}{\|\mathbf a\|^3}.
\]

参考曲率：

\[
\kappa_0=\frac1R.
\]

弯曲应变：

\[
\chi=\kappa-\kappa_0.
\]

弯曲内力：

\[
\mathbf F_b^{\mathrm{int}}
=\int\mathbf B_b^T EI\chi\,dS.
\]

完成后，总内力为：

\[
\mathbf F^{\mathrm{int}}
=\mathbf F_m^{\mathrm{int}}
+\mathbf F_b^{\mathrm{int}}.
\]

### 通过标准

- `d=0` 时膜内力和弯曲内力都接近 0；
- 刚体平移不产生内力；
- 刚体有限转动不应产生虚假应变能；
- 数值切线近似对称；
- 小载荷结果退化到线性解。

---

## 阶段 7：推导并实现解析切线

最后才把中心差分切线替换为：

\[
\mathbf K_t
=\mathbf K^{\mathrm{mat}}
+\mathbf K^{\mathrm{geo}}.
\]

实现顺序：

1. 先推导 membrane 的 `Kmat`；
2. 再推导 membrane 的 `Kgeo`；
3. 与完整内力的中心差分 Jacobian 对比；
4. 再推导 bending 切线；
5. 最后关闭数值切线。

### 通过标准

解析切线与中心差分切线满足：

`norm(Kanalytic-Kfd,'fro')/max(norm(Kfd,'fro'),eps) < tolerance`

其中初始学习阶段可以先取 `tolerance = 1e-5`。

---

## 阶段 8：最终验证与代码整理

需要完成的检查：

1. 零位移、刚体平移和刚体转动测试；
2. 小载荷线性退化测试；
3. 载荷步收敛性测试；
4. 节点数和支撑域敏感性测试；
5. 数值/解析切线对比；
6. 位移、应变、曲率和残量输出；
7. 与完整 nonlinear 参考程序对照；
8. 最后再进入 Raknes 4.2 端弯矩算例。

## 每个阶段的工作方式

每次只推进一个阶段：

1. 先指出原线性程序中需要修改的具体代码块；
2. 先解释公式和变量尺寸；
3. 由学习者完成代码；
4. 检查代码和运行结果；
5. 记录通过标准；
6. 再进入下一阶段。

第一项实际任务是阶段 0：运行线性基准并记录位移、相对误差、刚度对称性和约束残量。
