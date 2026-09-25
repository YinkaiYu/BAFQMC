# 模型、晶格与物理约定

BAFQMC 在有限温度虚时路径积分中对连续 Hubbard–Stratonovich（HS）场进行采样。对于每个场构型，通过二次单粒子传播对玻色子求迹。HS 解耦后的二次型问题的时间反演对称性（TRS）或反射正性（RP）确保权重非负，包括在阻挫跃迁情形下。所提供的实现处理粒子数守恒玻色子和格点配对两种情况。

## 求解器实现的哈密顿量

我们采用论文的算符与符号约定：两个玻色子组分 $`\hat b_i`$ 和
$`\hat c_i`$，其中 $`\hat n_{b,i}=\hat b_i^+\hat b_i`$ 和
$`\hat n_{c,i}=\hat c_i^+\hat c_i`$。主要基准哈密顿量为

```math
\begin{aligned}
\hat H={}&\hat H_t+\hat H_U+\hat H_\Delta,\\
\hat H_t={}&t\sum_{\langle ij\rangle}
\left(\hat b_i^+\hat b_j+\hat b_j^+\hat b_i
+\hat c_i^+\hat c_j+\hat c_j^+\hat c_i\right),\\
\hat H_U={}&U\sum_i(\hat n_{b,i}-\hat n_{c,i})^2,\\
\hat H_\Delta={}&-\sum_i\left(\Delta\,\hat b_i^+\hat c_i^+
+\Delta^*\,\hat b_i\hat c_i\right).
\end{aligned}
```

系综为巨正则系综，逆温度为 $`\beta`$：

```math
\hat N=\sum_i(\hat n_{b,i}+\hat n_{c,i}),\qquad
\hat H_\mu=\hat H-\mu\hat N,\qquad
Z=\mathrm{Tr}e^{-\beta\hat H_\mu}.
```

粒子数守恒求解器取 $`\Delta=0`$。配对求解器接受实数 $`\Delta`$，并保留正常和反常 Green 函数。其内部算符相位为 $`\hat b_{\mathrm{code}}=\hat b_{\mathrm{paper}}`$，
$`\hat c_{\mathrm{code}}=-\hat c_{\mathrm{paper}}`$，因此在相同输入 $`\Delta`$ 下，代码中的配对项带正号。这是相同的物理模型。
实现中的 Nambu 顺序为 $`(b,c,b^+,c^+)`$（在该代码基底下）。

两个求解器还实现了补充材料中使用的两个密度通道，将 $`\hat H_U`$ 替换为

```math
\begin{aligned}
\hat H_U={}&\sum_i\left[
U_1(\hat n_{b,i}-\hat n_{c,i})^2
+U_2(\hat n_{b,i}+\hat n_{c,i})^2\right]\\
={}&\sum_i\left[(U_1+U_2)(\hat n_{b,i}^2+\hat n_{c,i}^2)
+2(U_2-U_1)\hat n_{b,i}\hat n_{c,i}\right].
\end{aligned}
```

主模型设 $`U_1=U`$，$`U_2=0`$。补充材料固定 $`U_1=1`$ 并扫描 $`U_2\leq0`$。

平方密度包含其线性粒子数项。若从以 $`n(n-1)`$ 形式写出的模型转换，需将由此产生的化学势偏移同时代入 BAFQMC 和 ED。

## 三角晶格

当前几何结构为每个元胞一个格点，带有周期边界条件：

```math
\begin{gathered}
\mathbf a_1=(1,0),\qquad
\mathbf a_2=\left(\frac12,\frac{\sqrt3}{2}\right),\\
\mathbf r_{x,y}=x\mathbf a_1+y\mathbf a_2,\qquad
N_s=L_xL_y,\\
(x,y)\equiv(x+L_x,y)\equiv(x,y+L_y).
\end{gathered}
```

每个元胞提供三条正向键

```math
(x,y)\longrightarrow(x+1,y),\quad
(x,y)\longrightarrow(x,y+1),\quad
(x,y)\longrightarrow(x-1,y+1).
```

其厄米共轭提供反向。在非常小的周期元胞上，不同键可能连接同一对格点或回到原格点；它们的贡献被叠加。独立的 ED 和固定场测试明确使用了这种键的多重性。

在已实现的正跃迁 $`t=1`$ 下，三角环是阻挫的。色散关系为

```math
\varepsilon(\mathbf k)=2t\left[
\cos(\mathbf k\cdot\mathbf a_1)
+\cos(\mathbf k\cdot\mathbf a_2)
+\cos\bigl(\mathbf k\cdot(\mathbf a_2-\mathbf a_1)\bigr)\right].
```

倒格矢基矢和一个能带极小值为

```math
\begin{gathered}
\mathbf b_1=2\pi\left(1,-\frac1{\sqrt3}\right),\qquad
\mathbf b_2=2\pi\left(0,\frac2{\sqrt3}\right),\\
\mathbf b_i\cdot\mathbf a_j=2\pi\delta_{ij},\qquad
\mathbf K=\frac23\mathbf b_1+\frac13\mathbf b_2
=\left(\frac{4\pi}{3},0\right).
\end{gathered}
```

允许的动量为 $`\mathbf k=(m/L_x)\mathbf b_1+(n/L_y)\mathbf b_2`$。
现有的 $`K`$ 点估计量要求两个长度均为 3 的倍数。在其他尺寸下，这些通道的 Fortran 输出为零；若进行新的研究，请选择相容的动量并实现其估计量。密度和能量测量在其他晶格尺寸上仍然可用。

## 参数与当前接口

| 物理量 | 现有控制方式 |
| --- | --- |
| 晶格尺寸 $`L_x,L_y`$ | 运行时输入 `Nlx,Nly`；BAFQMC 接受不同的正整数晶格长度 |
| 温度 $`\beta`$ 和时间步长 $`\Delta\tau=\beta/L_\tau`$ | 运行时输入 `Beta,Ltrot` |
| 密度相互作用 $`U_1,U_2`$ 和化学势 $`\mu`$ | 运行时输入 `U1,U2,mu` |
| 实数格点配对 $`\Delta`$ | 配对求解器中的运行时输入 `RDelta` |
| 跃迁强度 $`t`$ | 在每个求解器的 `src/calc_basic.f90` 中由 `Params_set` 设为 `RT=1.d0`；若需改变，请修改实现并重新编译 |
| 晶体几何、键振幅、额外相互作用或配对模式 | 模型开发；使用[扩展指南](./model-development.md) |

仅有 JSON 的 `t` 字段无法改变当前 Fortran 的跃迁强度。所提供的粒子数守恒相互作用 ED 工作流针对 $`3\times3`$ 基准晶格；配对 ED 实现在其基底截断范围内接受小的 $`L_x\times L_y`$ 元胞。将 BAFQMC 计算扩展到新几何结构时，需选择并实现相应的参考计算。

## 从相互作用到辅助场

令 $`n_+=n_b+n_c`$，$`n_-=n_b-n_c`$。局域 HS 恒等式为

```math
\begin{aligned}
e^{-\Delta\tau U_1 n_-^2}
&=\int\frac{d\phi_1}{\sqrt{2\pi}}e^{-\phi_1^2/2}
  e^{\mathrm i\sqrt{2U_1\Delta\tau}\,\phi_1n_-},\\
e^{-\Delta\tau U_2 n_+^2}
&=\int\frac{d\phi_2}{\sqrt{2\pi}}e^{-\phi_2^2/2}
  e^{\sqrt{-2U_2\Delta\tau}\,\phi_2n_+}.
\end{aligned}
```

对于粒子数守恒传播，两个组分对应共轭的在位势。因此 $`B_{c,\ell}=\overline{B_{b,\ell}}`$，代码可通过对 $`b`$ 扇区取共轭得到 $`c`$ 扇区的 Green 函数。配对实现对完整的 Nambu 高斯迹求值，包括正规序标量因子和行列式平方根权重。局域场更新、稳定化传播和 Wick 估计量在 Fortran 中实现；Python 负责 ED、计算方案和分析。

对于主模型 $`U_2=0`$，充分条件 $`\mu<-3t-|\Delta|`$ 在整个辅助场域内确保迹有限；所有主要基准点均满足该条件。补充材料中具有吸引相互作用 $`U_2<0`$ 的 benchmark 使用其单独指定的有限占据比较参考。对于新的哈密顿量，请按[模型开发指南](./model-development.md)中的说明，建立迹域及相应的 HS 对称性与实现。

## 基准模型与可观测量

| 扫描 | 哈密顿量参数 | 几何结构与系综 |
| --- | --- | --- |
| 主相互作用扫描 | $`U_1=U`$，$`U_2=0`$，$`\Delta=0`$ | $`3\times3`$，$`t=1`$，$`\beta=4`$，$`\mu=-3.5`$ |
| 主配对扫描 | $`U_1=1`$，$`U_2=0`$，变化 $`\Delta`$ | $`3\times3`$，$`t=1`$，$`\beta=4`$，$`\mu=-5`$ |
| 补充密度通道扫描 | $`U_1=1`$，变化 $`U_2\leq0`$，$`\Delta=0`$ | $`3\times3`$，$`t=1`$，$`\beta=1`$，$`\mu=-7`$ |

论文中配对项带负号。论文算符与实现算符的关系为
$`b_{\mathrm{code}}=b_{\mathrm{paper}}`$ 和
$`c_{\mathrm{code}}=-c_{\mathrm{paper}}`$。密度、能量以及两个基准结构因子在此相位变换下不变；反常配对振幅的符号改变。

对于 $`\hat n_i=\hat n_{b,i}+\hat n_{c,i}`$，主要基准可观测量为

```math
\begin{aligned}
\rho&=\frac{\langle\hat N\rangle}{N_s},\qquad E=\langle\hat H\rangle,\\
S_{\mathrm{SF}}(\mathbf k)&=\frac1{N_s^2}\sum_{ij}
e^{\mathrm i\mathbf k\cdot(\mathbf r_i-\mathbf r_j)}
\left\langle\hat b_i^+\hat b_j+\hat c_i^+\hat c_j\right\rangle,\\
S_{\mathrm{DW}}(\mathbf k)&=\frac1{N_s^2}\sum_{ij}
e^{\mathrm i\mathbf k\cdot(\mathbf r_i-\mathbf r_j)}
\left\langle\hat n_i\hat n_j\right\rangle.
\end{aligned}
```

[输出可观测量指南](./observables.md)定义了每个当前物理输出量及其归一化方式。`energy_density` 为 $`E/N_s`$，不含化学势项；论文绘制的是 $`-E=-N_s\,\texttt{energy\_density}`$。代码的反常振幅为
$`\texttt{pair\_equal}=N_s^{-1}\sum_i\langle b_i^+c_i^++b_ic_i\rangle_{\mathrm{code}}`$，
因此每格点配对能为 $`\Delta\,\texttt{pair\_equal}`$。
论文的反常振幅为该值的相反数。

工作流报告均值、基于分块的均值标准误（SEM）、参考值及差值。主配对 ED 设置 `nmax=3,ncut=4` 包含 7297 个基底态；粒子数守恒 ED 累积粒子数和平移扇区。非相互作用参考值也直接从自由玻色色散计算得到。每次比较均记录其采样参数和 ED 截断值。

参见[基准输入与设置](./benchmarks-paper-readme.md)、[科学测试](./testing.md)，以及详细的[粒子数守恒](https://github.com/YinkaiYu/BAFQMC/blob/main/docs/solvers/number_conserving/physics.md)和[配对](https://github.com/YinkaiYu/BAFQMC/blob/main/docs/solvers/pairing/physics.md)实现指南。
