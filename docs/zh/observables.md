# 物理观测量与输出文件

两个求解器都在周期边界单轨道三角晶格上测量以下物理量。
定义采用论文的算符约定，$`N_s=L_xL_y`$（代码中为 `Lq`）：

```math
\hat n_{b,i}=\hat b_i^+\hat b_i,\qquad
\hat n_{c,i}=\hat c_i^+\hat c_i,\qquad
\hat N_b=\sum_i\hat n_{b,i},\qquad
\hat N_c=\sum_i\hat n_{c,i}.
```

期望值均为巨正则热平均。文件名中的 `up` 和 `do` 分别对应 $`b`$ 和 $`c`$。
除特别标注 **仅配对求解器** 外，以下所有文件均由两个求解器写出。
配对求解器的反常振幅有特殊的相位约定，详见下文；密度、跃迁和三个结构因子
与论文的约定直接一致。

## 密度与粒子数矩

定义平均密度和格点二阶矩求和为

```math
\rho_i=\langle\hat n_{b,i}+\hat n_{c,i}\rangle,\qquad
\rho=\frac{\langle\hat N_b+\hat N_c\rangle}{N_s},\qquad
M_b=\sum_i\langle\hat n_{b,i}^2\rangle,\quad
M_c=\sum_i\langle\hat n_{c,i}^2\rangle.
```

| 输出文件 | 定义 |
| --- | --- |
| `density_up`, `density_do` | $`\langle\hat N_b\rangle/N_s`$，$`\langle\hat N_c\rangle/N_s`$ |
| `density_total` | $`\rho`$ |
| `density`（**仅配对求解器**） | $`\rho`$ 的别名 |
| `density_site_total` | $`N_s`$ 个格点密度 $`\rho_i`$ |
| `num_up`, `num_do` | $`\langle\hat N_b\rangle`$，$`\langle\hat N_c\rangle`$ |
| `numsquare_up`, `numsquare_do` | $`\langle\hat N_b^2\rangle`$，$`\langle\hat N_c^2\rangle`$ |
| `onsite_n2_up`, `onsite_n2_do` | $`M_b`$，$`M_c`$ |
| `local_numsquare`（**仅配对求解器**） | $`(M_b+M_c)/N_s`$ |
| `doubleOcc` | $`D=N_s^{-1}\sum_i\langle\hat n_{b,i}\hat n_{c,i}\rangle`$ |
| `squareOcc` | $`Q=(2N_s)^{-1}\sum_i\langle\hat n_{b,i}(\hat n_{b,i}-1)+\hat n_{c,i}(\hat n_{c,i}-1)\rangle`$（当前模型） |

`numsquare_*` 包含不同格点间的关联；`onsite_n2_*` 只含格点自身项。
这些文件存储总量，而 `doubleOcc`、`squareOcc` 和 `local_numsquare` 为每格点值。

辅助场估计量 `squareOcc` 的精确表达式为

```math
Q_\phi=\frac1{N_s}\sum_i\mathrm{Re}\left[
\langle\hat n_{b,i}\rangle_\phi^2+
\langle\hat n_{c,i}\rangle_\phi^2\right].
```

对于当前的哈密顿量，$`\hat N_b-\hat N_c`$ 的守恒使同组分反常缩并为零，
因此 Wick 定理给出表中的正规序定义。含同组分配对的模型需要额外的反常缩并。
配对求解器的 `local_numsquare` 和 `onsite_n2_*` 估计量已显式包含这些缩并项。

密度分布的逆参与比（IPR）是后处理量，求解器无单独输出文件：

```math
\mathrm{IPR}_\rho=
\frac{\sum_i\rho_i^2}{\left(\sum_i\rho_i\right)^2}.
```

它使用平均后的密度分布，对于均匀非零分布等于 $`1/N_s`$，
与粒子数二阶矩不同。

## 能量与反常配对振幅

物理哈密顿量 $`\hat H=\hat H_t+\hat H_U+\hat H_\Delta`$，其中

```math
\hat H_t=t\sum_{\langle ij\rangle}\left(
\hat b_i^+\hat b_j+\hat b_j^+\hat b_i+
\hat c_i^+\hat c_j+\hat c_j^+\hat c_i\right),
```

```math
\hat H_U=U_1\sum_i(\hat n_{b,i}+\hat n_{c,i})^2+
U_2\sum_i(\hat n_{b,i}-\hat n_{c,i})^2,
```

```math
\hat H_\Delta=-\Delta\sum_i
\left(\hat b_i^+\hat c_i^++\hat b_i\hat c_i\right).
```

配对求解器接受实数 $`\Delta`$；粒子数守恒求解器取 $`\Delta=0`$。
论文主模型取 $`U_1=0`$，$`U_2=U`$。
$`\hat H_t`$ 中每条最近邻键只计一次，显式写出两个跃迁方向。
定义 $`e_t=\langle\hat H_t\rangle/N_s`$，
$`e_U=\langle\hat H_U\rangle/N_s`$，
$`e_\Delta=\langle\hat H_\Delta\rangle/N_s`$。

配对求解器使用

```math
\hat b_{i,\mathrm{code}}=\hat b_i,\qquad
\hat c_{i,\mathrm{code}}=-\hat c_i,
```

因此其二次型矩阵中配对项符号为正。由此

```math
\begin{aligned}
P_{\mathrm{paper}}&=\frac1{N_s}\sum_i
\langle\hat b_i\hat c_i+\hat b_i^+\hat c_i^+\rangle,\\
P_{\mathrm{code}}&=\frac1{N_s}\sum_i
\langle\hat b_{i,\mathrm{code}}\hat c_{i,\mathrm{code}}+
\hat b_{i,\mathrm{code}}^+\hat c_{i,\mathrm{code}}^+\rangle
=-P_{\mathrm{paper}},\\
e_\Delta&=\Delta P_{\mathrm{code}}=-\Delta P_{\mathrm{paper}}.
\end{aligned}
```

| 输出文件 | 定义 |
| --- | --- |
| `kinetic` | $`e_t`$，含 $`t`$ 和两种组分 |
| `interaction_energy_density` | $`e_U=(U_1+U_2)(M_b+M_c)/N_s+2(U_1-U_2)D`$ |
| `pair_equal`（**仅配对求解器**） | $`P_{\mathrm{code}}`$；取负号得论文的配对振幅 |
| `pairing_energy_density`（**仅配对求解器**） | $`e_\Delta=\Delta P_{\mathrm{code}}`$ |
| `energy_density` | $`e=\langle\hat H\rangle/N_s=e_t+e_U+e_\Delta`$ |
| `chemical_energy_density`（**仅配对求解器**） | $`e_\mu=-\mu\rho`$ |
| `grand_energy_density`（**仅配对求解器**） | $`\langle\hat H-\mu(\hat N_b+\hat N_c)\rangle/N_s=e-\mu\rho`$ |

`energy_density` 不含化学势项。总物理能量 $`E=N_s e`$；论文绘制 $`-E`$。
乘以 $`N_s`$ 同时也将标准误差乘以 $`N_s`$。
粒子数守恒求解器不单独写出化学势贡献或巨正则能量文件；它们由相同的恒等式推导得到。

## K 点和 Γ 点结构因子

三角晶格原胞基矢和测量动量为

```math
\mathbf a_1=(1,0),\qquad
\mathbf a_2=\left(\frac12,\frac{\sqrt3}{2}\right),\qquad
\mathbf K=\left(\frac{4\pi}{3},0\right),\qquad
\boldsymbol\Gamma=(0,0).
```

三个结构因子均采用 $`N_s^{-2}`$ 归一化：

```math
S_{\mathrm{SF}}(\mathbf K)=\frac1{N_s^2}\sum_{ij}
e^{i\mathbf K\cdot(\mathbf r_i-\mathbf r_j)}
\langle\hat b_i^+\hat b_j+\hat c_i^+\hat c_j\rangle,
```

```math
S_{\mathrm{DW}}(\mathbf K)=\frac1{N_s^2}\sum_{ij}
e^{i\mathbf K\cdot(\mathbf r_i-\mathbf r_j)}
\langle(\hat n_{b,i}+\hat n_{c,i})(\hat n_{b,j}+\hat n_{c,j})\rangle,
```

```math
S_{\mathrm{PSF}}(\boldsymbol\Gamma)=\frac1{N_s^2}\sum_{ij}
\langle\hat b_i^+\hat c_i^+\hat c_j\hat b_j\rangle.
```

| 输出文件 | 可观测量 | 可用条件 |
| --- | --- | --- |
| `sf_K` | $`S_{\mathrm{SF}}(\mathbf K)`$ | $`L_x`$ 和 $`L_y`$ 均为 3 的倍数 |
| `dw_K` | $`S_{\mathrm{DW}}(\mathbf K)`$ | $`L_x`$ 和 $`L_y`$ 均为 3 的倍数 |
| `psf_Gamma` | $`S_{\mathrm{PSF}}(\boldsymbol\Gamma)`$ | 所有支持的晶格尺寸 |

`dw_K` 测量总密度。这些是完整的关联函数，包含非连通贡献，
不减去平均值的乘积。配对结构因子在代码到论文的相位旋转下不变。

当晶格无法表示 K 点时，可执行程序仍会写出 `sf_K` 和 `dw_K`，
但填充零占位符。对该晶格应将这些通道视为不可用，不影响 `psf_Gamma`。

## 密度关联文件族

写入程序对平移平均的密度关联做正号傅里叶变换

```math
C_{bb}(\mathbf q)=\frac1{N_s^2}\sum_{ij}
e^{i\mathbf q\cdot(\mathbf r_i-\mathbf r_j)}
\langle\hat n_{b,i}\hat n_{b,j}\rangle,
```

类似地定义 $`C_{cc}`$ 和 $`C_{bc}`$（替换相应的密度算符）。
平移平均贡献 $`1/N_s`$，傅里叶变换再贡献 $`1/N_s`$。
当前写入程序只选取 $`\mathbf q=\boldsymbol\Gamma`$：

| 输出文件 | 写入量 |
| --- | --- |
| `den_upup_sub11` | $`C_{bb}(\boldsymbol\Gamma)=\langle\hat N_b^2\rangle/N_s^2`$ |
| `den_dodo_sub11` | $`C_{cc}(\boldsymbol\Gamma)=\langle\hat N_c^2\rangle/N_s^2`$ |
| `den_updo` | $`C_{bc}(\boldsymbol\Gamma)=\langle\hat N_b\hat N_c\rangle/N_s^2`$ |

`sub11` 后缀表示单轨道的两个关联指标。这些文件每个 bin 包含一个复数值，
不含动量坐标或完整动量网格，也不减去非连通贡献。
特别地，`den_updo` 包含不同格点间的关联，而 `doubleOcc` 只含格点自身关联。

## 行、列与不确定度

文件在运行的工作目录中追加写入。每次独立运行请使用新目录。
每个等时行是一个 bin 均值：局域扫描对虚时片和蒙特卡洛扫描的测量取平均，
然后由 0 号进程写出 MPI 进程的平均值。
对于 `Nsweep` 次双向扫描，每个进程每个 bin 贡献 `2 * Nsweep * Ltrot` 个等时测量。

| 文件类型 | 每行的列 |
| --- | --- |
| 标量密度、矩、能量和 `pair_equal` | 一个实数值 |
| `density_site_total` | $`N_s`$ 个实数值，每格点一个 |
| `sf_K`、`dw_K`、`psf_Gamma` 和 `den_*` | 两个值：实部、虚部 |

对于 `density_site_total`，零基格点指标 $`i=x+L_x y`$ 对格点排序，
$`x=0,\ldots,L_x-1`$ 变化最快。原始文件不存储列标签。
上述精确结构因子和密度关联在理论上为实数；虚部列保留有限样本估计量，
作为一致性诊断使用。

后处理对保留的 bin 序列进行分块以估计均值的不确定度。
对于 $`B`$ 个分块均值 $`\bar x_j`$：

```math
\bar x=\frac1B\sum_{j=1}^B\bar x_j,\qquad
\mathrm{SEM}(\bar x)=
\sqrt{\frac{\sum_{j=1}^B(\bar x_j-\bar x)^2}{B(B-1)}}.
```

报告的 `stderr` 是该标准误差，而非原始 bin 的标准差。复数列分别分析。
原始样本在本地生成；小型处理后的表格保留均值、均值标准误及用于重建统计量的分块均值。
论文工作流详见[复现指南](./benchmarks-paper-readme.md)。

## 已启用测量与诊断

等时测量从第一个输出 bin 开始。`Nthermal` 控制非等时路径的进入时机；
它不丢弃等时行。辅助场预热和任何分析时的 bin 丢弃是独立步骤。

两个求解器当前的 `obser_tau.f90` 估计量和 `m_write_obs_tau` 写入程序均为空。
因此开启 `is_tau` 不产生任何非等时的物理可观测量、磁化率或松原频率数据。
配对求解器会累积内部的单粒子关联，但其 `green` 写入调用已被注释掉。
通用傅里叶输出辅助方法本身并不开启额外的输出文件。

粒子数守恒求解器的文件 `pole_z`、`pole_distance`、`pole_x`、
`green_spectral_radius`、`green_smax` 和 `log_weight` 描述构型权重和数值条件数。
与等时可观测量行不同，这些文件**每个 bin 每个 MPI 进程**包含一个构型样本，
按进程顺序排列。其定义参见[粒子数守恒物理指南](https://github.com/YinkaiYu/BAFQMC/blob/main/docs/solvers/number_conserving/physics.md#continuous-pole-diagnostics)。
`info.txt` 中的接受率、相位和矩阵稳定性统计量也是诊断量。
日志、种子和场构型记录运行时状态。

## 实现地图

| 职责 | 粒子数守恒求解器 | 配对求解器 |
| --- | --- | --- |
| Wick 估计量和归一化 | [obser_equal.f90](https://github.com/YinkaiYu/BAFQMC/blob/main/src/number_conserving/src/obser_equal.f90) | [obser_equal.f90](https://github.com/YinkaiYu/BAFQMC/blob/main/src/pairing/src/obser_equal.f90) |
| 已启用文件、列布局、MPI 平均 | [fourier_trans.f90](https://github.com/YinkaiYu/BAFQMC/blob/main/src/number_conserving/src/fourier_trans.f90) | [fourier_trans.f90](https://github.com/YinkaiYu/BAFQMC/blob/main/src/pairing/src/fourier_trans.f90) |
| 晶格排序和傅里叶归一化 | [lattice.f90](https://github.com/YinkaiYu/BAFQMC/blob/main/src/number_conserving/src/lattice.f90) | [lattice.f90](https://github.com/YinkaiYu/BAFQMC/blob/main/src/pairing/src/lattice.f90) |
| Bin 循环和测量标志 | [main.f90](https://github.com/YinkaiYu/BAFQMC/blob/main/src/number_conserving/src/main.f90) | [main.f90](https://github.com/YinkaiYu/BAFQMC/blob/main/src/pairing/src/main.f90) |
| 非等时估计量状态 | [obser_tau.f90](https://github.com/YinkaiYu/BAFQMC/blob/main/src/number_conserving/src/obser_tau.f90) | [obser_tau.f90](https://github.com/YinkaiYu/BAFQMC/blob/main/src/pairing/src/obser_tau.f90) |

扩展时，请同步更新每个估计量、其归一化以及已启用的写入程序。
[模型扩展指南](./model-development.md)提供了对应哈密顿量、基组和更新的修改地图。
