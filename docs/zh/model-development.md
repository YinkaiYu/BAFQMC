# 实现新模型

BAFQMC 可以在本仓库所发布的两个三角晶格模型之外进一步开发。本指南将新的哈密顿量映射到实际的求解器、ED 和测量代码。内容涵盖晶格、跃迁、配对、组分和相互作用的修改。若要对已实现的模型进行参数扫描，请使用[计算方案](./agent-workflows.md)。[可观测量参考](./observables.md)定义了当前测量的算符和输出归一化，扩展工作可在此基础上进行。

从研究人员请求的哈密顿量出发，将扩展工作贯穿至一个可运行的示例和独立的数值检验。下面给出的各向异性晶格示例是一个**实现蓝图**；当前可执行文件没有通用晶格或跃迁矩阵的输入标志。

## 指定格点、算符与哈密顿量

区分元胞、物理格点、组分模式和 Nambu 分量。对于每个元胞含 $`n_{\mathrm{sub}}`$ 个格点的晶体，

```math
N_{\mathrm{cell}}=L_xL_y,\qquad
N_s=n_{\mathrm{sub}}N_{\mathrm{cell}},\qquad
\mathbf r_{\mathbf R a}=R_1\mathbf a_1+R_2\mathbf a_2+\boldsymbol\delta_a.
```

其中 $`a`$ 是子晶格指标，$`\boldsymbol\delta_a`$ 是其在元胞内的位置。有 $`F`$ 个玻色组分时，共有 $`M=F N_s`$ 个湮灭算符。完整的 Nambu 表示含 $`2M`$ 个分量。现有的粒子数守恒求解器传播一个组分并重建其共轭组分；因此对于当前模型，其矩阵维数为 $`N_s`$。

写出格点/组分指标约定和明确的键列表，包括周期镜像位移。说明每条键是以厄米共轭形式出现一次，还是以两条有向条目的形式出现。小周期团簇可能有重复的镜像键和自镜像键：对其振幅求和与对端点去重是不同的哈密顿量。当前代码对镜像贡献求和。

一个有用的通用二次型规范——采用论文的负号配对约定——为

```math
H_0
=\sum_{\alpha\beta} a_\alpha^+ h^{(0)}_{\alpha\beta}a_\beta
-\frac12\sum_{\alpha\beta}
\left(\Delta_{\alpha\beta}a_\alpha^+a_\beta^+
+\Delta_{\alpha\beta}^{*}a_\beta a_\alpha\right)+C,
\qquad h^{(0)}=\bigl(h^{(0)}\bigr)^\dagger,\quad \Delta=\Delta^T.
```

物理哈密顿量及其巨正则对应量为

```math
H=H_0+H_{\mathrm{int}},\qquad
H_\mu=H-\mu N,\qquad N=\sum_\alpha a_\alpha^+a_\alpha,
\qquad h=h^{(0)}-\mu I.
```

联合指标 $`\alpha,\beta`$ 包含物理格点和组分，$`h`$ 是用于热迹的粒子数守恒二次项矩阵。厄米跃迁要求反向振幅为复共轭。玻色配对在两个联合指标交换下对称；因子 $`1/2`$ 避免重复计数。明确说明物理配对相位。当前求解器使用正实数的组分间配对系数，与论文的负号系数之间的关系为 $`c_{\mathrm{code}}=-c_{\mathrm{paper}}`$。

用算符形式写出相互作用，包括其线性项和常数项。例如，$`n^2=n(n-1)+n`$，因此用一个替换另一个时也会改变化学势的贡献。在整个实现过程中，分别定义物理能量和巨正则能量。对于论文主模型的扩展，用 $`b,c`$ 表示两个组分，用 $`U=U_2`$、$`U_1=0`$ 表示相对密度相互作用。

## 推导解耦方案及其符号保护

对于写成厄米二次通道之和的相互作用 $`H_{\mathrm{int}}=\sum_\ell g_\ell Q_\ell^2`$，推导每个时间步长为 $`\Delta\tau`$ 的 Trotter 因子中所使用的 HS 表示。连续高斯通道具有恒等式

```math
e^{-\Delta\tau g Q^2}
=\int_{-\infty}^{\infty}\frac{d\phi}{\sqrt{2\pi}}
e^{-\phi^2/2}
e^{\sqrt{-2\Delta\tau g}\,\phi Q}.
```

当 $`g<0`$ 时，场系数为实数；当 $`g>0`$ 时，场系数为虚数。若各通道或动能因子不对易，则记录其顺序和 Trotter 近似。将所有正规序移位和标量因子纳入权重，而非将全部信息放入矩阵传播子中。

对于每个独立的场构型，在解耦**之后**建立适用的时间反演对称性（TRS）、反射正性（RP）或共轭扇区构造。写出对称算符或反射、其在格点/组分/Nambu 指标上的作用，以及对每个生成的二次因子的条件。对于 TRS 路线，记录反幺正算符的平方以及所用权重非负判据的所有附加假设。对于 RP 路线，指定反射所交换的两个子系统，以及跨越该划分的耦合所需的形式和符号。原始相互作用本身的对称性并不能检验这些构型级别的条件。

最简单的现有构造给出了一个明确的例子。对于两个组分的等实跃迁，$`U_1\leq0`$ 和 $`U_2\geq0`$ 在**相同**场值下生成共轭的单组分因子。若 $`B_c(\phi)=B_b(\phi)^*`$，则在收敛热迹域内，

```math
w(\phi)=p_{\mathrm{HS}}(\phi)
\left|\det\left[I-B_b(\phi)\right]\right|^{-2}\geq0.
```

粒子数守恒代码直接使用此关系。具有不同组分跃迁、组分混合或新 HS 通道的模型，在该关系改变时，需要自己的矩阵表示和权重推导。对新行列式取绝对值并不能证明其符号保护。对于复数跃迁，该快捷方式要求共轭的粒子数守恒块 $`h_c=h_b^*`$ 以及共轭 HS 因子。两个组分具有相同通量的等复数跃迁通常不满足该条件。将一个组分的通量替换为其相反数会改变所请求的哈密顿量；请使用适用的对称性构造来分析实际的同通量模型。若所提议的哈密顿量不属于已建立的类，则指明哪个条件受到影响，并相应地开发所请求的模型；不要悄然改变其物理耦合以恢复旧的类。

若所请求的模型需要复数权重采样，则需显式实现该估计量。配对求解器的局域比率-相位诊断不是相位重加权：其当前测量是普通平均值。复数权重计算需要完整的构型相位 $`w(\phi)=|w(\phi)|e^{i\theta(\phi)}`$ 及相应的估计量

```math
\langle O\rangle
=\frac{\langle O(\phi)e^{i\theta(\phi)}\rangle_{|w|}}
{\langle e^{i\theta(\phi)}\rangle_{|w|}}.
```

这需要包含相位重加权的测量和分析，包括 $`\theta`$ 的标量和分支贡献，或对同一物理哈密顿量单独建立权重非负的构造。

## 建立物理迹与矩阵表示

符号保护和迹收敛是构造的独立部分。指定巨正则物理模型具有有限热迹的参数区域，并建立解耦后二次迹所需的域。有限 ED 占据截断定义了一个有限维参考；向未截断模型的收敛是进一步的计算。

对于二次模型，厄米玻色能量矩阵的正定性是有限热迹的一个有用充分条件。上述二次形式的一个充分下界为

```math
\lambda_{\min}(h)>\|\Delta\|_2.
```

仅凭实数 Bogoliubov 频率并不能建立该能量矩阵的正定性。对于当前具有均匀在位配对和实数跃迁的 $`U_1=0`$ 构造，相应的充分界为

```math
\mu<\varepsilon_{\min}-|\Delta|,
```

其中 $`\varepsilon_{\min}`$ 是跃迁能带的最小值。相对密度 HS 因子是幺正的。三角晶格 benchmark 有 $`\varepsilon_{\min}=-3t`$；不同的图需要自己的能带界。新的非幺正 HS 通道需要相应的迹域分析。

保持二次能量矩阵与被取指数的玻色对易子矩阵的区别。在完整 Nambu 形式中，即使对于厄米物理哈密顿量，后者通常也是非厄米的。配对求解器的 `exp_general_matrix` 使用通用复数对角化。新的配对矩阵必须保留正确的粒子/空穴符号和共轭；对该对易子矩阵使用厄米本征求解器会改变传播子。

配对表示还有一个标量正规序贡献。对于当前的总密度 HS 场 $`x=\sqrt{-2U_1\Delta\tau}\,\phi_1`$，标量因子为 $`e^{-x}`$，一次提议对 `ratio_constant` 贡献 $`e^{-(x'-x)}`$。对不同的二次生成元或 Nambu 约定，需重新推导此因子。

局域行列式因子当前使用 `det_Pblock**(-0.5d0)`。代码记录相位并对模采样；它不维护显式的分支连续状态。为新配对模型建立与物理迹相关的平方根分支，以已知收敛极限为锚点，并检验包含标量的完整比率。若模型需要分支追踪，则需将该状态作为扩展的一部分加以实现和验证。仅改变行列式指数并不能在完整和简化的 Nambu 表示之间转换。

## 修改实际实现

下表指向当前启用的源文件。两个组件使用相似的名称，但矩阵维数和更新公式不同。

| 功能 | 文件与符号 | 所需协调 |
| --- | --- | --- |
| 输入与维数 | `src/<solver>/src/calc_basic.f90`：`read_input`、`Params_set` | 读取并广播新参数；区分元胞、格点、组分和扇区；将其记录在运行元数据中。`RT=1` 当前在此赋值。|
| 晶体与图 | `src/<solver>/src/lattice.f90`：`Lattice_make` | 更新 `L_bonds`、`LT_bonds`、格点/扇区指标映射、实空间/倒空间向量和 Fourier 相位。|
| 二次传播 | `src/<solver>/src/non_interact.f90`：`def_hamT`、`opT_set`；配对版 `exp_general_matrix` | 构建新的跃迁/配对矩阵及其逆时因子，保持物理约定和对易子约定的区别。|
| 场与通道调度 | `fields.f90`：`AuxConf_make`；`model.f90`：`Model_init`；`local_sweep.f90` | 更新场作用范围、通道分配、初始/重启 I/O 以及两个扫描方向。当前调度明确使用两个通道。|
| HS 因子 | `operator_Hubbard.f90`：`opU_set`、`opU_get_delta`、`opU_mmult_L/R` | 实现推导的算符、场测度、标量和局域矩阵变化。在配对代码中，耦合符号当前用于选择通道类型。|
| 比率与 Green 更新 | `localU.f90`：`LocalU_metro`；`multiply.f90`；`stabilization.f90` | 匹配更新作用范围和秩，推导行列式比率，并对照直接矩阵乘积验证稳定化传播。|
| 测量 | `obser_equal.f90`：`Obs_equal_calc`；`fourier_trans.f90` | 随哈密顿量更改动能、相互作用、配对和动量估计量；将定义传入文件输出和分析。|
| 计算与参考 | `src/<solver>/run_paper.py`、`src/<solver>/benchmarks/campaign*.py`、以下 ED 文件 | 生成一致的求解器/ED 输入，验证实现的模型，并报告每个新参数和参考截断。|

请直接通过[粒子数守恒组件](https://github.com/YinkaiYu/BAFQMC/blob/main/src/number_conserving/src/)和[配对组件](https://github.com/YinkaiYu/BAFQMC/blob/main/src/pairing/src/)查阅源码。特别注意：

- **子晶格数不是完整的晶格定义。** 当前的 `Norb=1` 键构造明确以目标轨道 1 为目标。Kagome 扩展需要三个物理子晶格、元胞内位置、元胞间键和一致的轨道指标。使用 `Lq`、`Ndim`、`Nsite`、`Norb`、`Nsub` 和 `Nbond` 审查分配和归一化。子晶格数、配位数和前向键存储数是不同的。`Nbond=3` 以及 `LT_bonds` 的空间/时间列当前固定为三角图。检查 Fourier 辅助函数（如 `m_write_k_3` 和 `m_write_reciprocal_3`）中的张量范围与其调用者是否一致；轨道关联维数必须遵循新的轨道索引。配对版的 `Nsec=4` 标记 $`(b,c,b^+,c^+)`$；它不是子晶格数。
- **二次项和相互作用的更改影响不同的更新作用范围。** 新的确定性跃迁或键配对可能保留现有的在位密度 HS 作用范围。键 HS 场、组分非对角通道或新相互作用会改变它。粒子数守恒的 `LocalU_metro` 应用秩一对角格点更新。配对版选择一个格点的四个扇区，并只使用其 $`4\times4`$ `Delta` 的对角项。仅向该数组添加非对角项并不能实现所需的 Woodbury 更新。
- **通道符号编码算符。** 配对版 `opU_set` 对负耦合选择总密度，对正耦合选择相对密度。相反符号或附加算符需要显式的通道定义；重命名 `Op_U1` 或更改 JSON 数字并不会改变该调度。
- **测量中明确包含旧哈密顿量。** 动能在 `L_bonds` 上对标量 `RT` 求和；配对能量使用在位的 `RDelta*pair_equal`。当传播发生变化时更新这些。当前配对实现提取完整的四扇区缩并；它不通过共轭 $`b`$ 来重构与之相互作用的 $`c`$ 扇区。
- **动量和归一化遵循新晶体。** 使用 $`e^{i\mathbf q\cdot(\mathbf r_{\mathbf R a}-\mathbf r_{\mathbf R'b})}`$（包含子晶格位置），并说明结果是子晶格指标上的矩阵还是求和后的物理结构因子。审查每个 $`N_s^{-1}`$ 和 $`N_s^{-2}`$ 因子。三角晶格 K 指标不是通用的有序态波矢；当前不可公度的 Fortran 运行会将其输出置零，而配对 ED 则忽略它。该零代表一个不可用的 K 估计量。定义并实现新模型的实际动量。

现有的配对含时可观测量例程 `Obs_tau_calc` 是一个占位符。涉及虚时关联的研究任务，除等时工作流外，还需要实现估计量及其传播检验。

## 扩展 ED 与分析以配合求解器

[粒子数守恒 ED 驱动](https://github.com/YinkaiYu/BAFQMC/blob/main/src/number_conserving/benchmarks/ed/EDtriangle_quspin_3x3.py)包含 `_hamiltonian_static`、`_hopping_two_species`、`_translation_permutations`、双组分基组构造和 K 点算符。其公开验证器当前将相互作用参考计算限制在 $`3\times3`$。引入其他尺寸或图时，请修改这些限制和对称块；不要保留新哈密顿量不再守恒的平移块。`src/number_conserving/benchmarks/campaign_analysis.py` 中的自由参考也包含三角色散，需要新的能谱。

对于配对，请更新 [geometry_triangle.py](https://github.com/YinkaiYu/BAFQMC/blob/main/src/pairing/benchmarks/ed/geometry_triangle.py)、通用驱动的 `build_basis`、`build_hamiltonian` 和可观测量算符，以及[张量 ED 实现](https://github.com/YinkaiYu/BAFQMC/blob/main/src/pairing/benchmarks/ed/ed_pairing_triangle_tensor.py)。当前基每个物理格点含两个组分。复数跃迁/配对还需要复数系数和合适的哈密顿量 dtype。保持物理格点数与 QuSpin 的 `basis.Ns`（即多体 Hilbert 空间维数）之间的区别。

当多个实现共存时，在新输入模式中添加模型标识符和明确的参数字段。扩展使用这些字段的读取器和验证器。新的清单标志只有在这些读取器和数值核心实现它之后才有效。保持生成的表格自描述：模型、几何、参数、单位、归一化、采样、SEM、ED 截断和能量约定。使用独立的示例和计算输出，以便原始三角形情形和新模型都可以运行。在定义或计算发生变化时有意更新参考数据，并在结果旁说明变化原因。

## 实现蓝图：各向异性三角跃迁

考虑一个具有相同两个组分和在位通道的请求扩展，但具有三个独立的实数跃迁振幅：

```math
H_{t}=\sum_{\mathbf R,\sigma}\sum_{\nu=1}^{3}t_\nu
\left(a_{\mathbf R\sigma}^+a_{\mathbf R+\mathbf d_\nu,\sigma}
+a_{\mathbf R+\mathbf d_\nu,\sigma}^+a_{\mathbf R\sigma}\right),
\qquad
\mathbf d_1=(1,0),\quad\mathbf d_2=(0,1),\quad\mathbf d_3=(-1,1).
```

位移对为整数元胞坐标。保留算法指南中的物理三角原胞基矢。色散变为

```math
\varepsilon(\mathbf k)=2\left[
t_1\cos(\mathbf k\cdot\mathbf a_1)
+t_2\cos(\mathbf k\cdot\mathbf a_2)
+t_3\cos\bigl(\mathbf k\cdot(\mathbf a_2-\mathbf a_1)\bigr)
\right].
```

对于等实数跃迁，该变化保留共轭组分关系和原始在位 HS 通道。对于 $`U_1=0`$ 和均匀实数在位配对，一个简单的充分热力学界为

```math
\mu<-2\left(|t_1|+|t_2|+|t_3|\right)-|\Delta|.
```

它由跃迁谱的下界推出，不必等于准确的能带最小值。作为具体的实现测试，使用 $`(t_1,t_2,t_3)=(1,0.8,0.35)`$、$`\mu=-5`$、$`\beta=1`$，首先令 $`U_1=U_2=0`$。测试 $`\Delta=0`$，以及对配对扩展测试 $`\Delta=0.1`$。然后在 $`U_2=0.7`$ 处检验现有的相对密度相互作用。这些是提议的扩展测试，不是捆绑的可运行情形。

实现此请求的 agent 应：

1. 在输入表示中添加并广播三个振幅；保留各向同性值 $`(1,1,1)`$ 作为原始模型的默认值。在 `def_hamT` 中将每个 `L_bonds(:,nu)` 条目连接到其振幅。
2. 使动能估计量使用相同的振幅和厄米共轭。在检验新二次因子的共轭扇区关系后，保持在位 HS 代数及其局域作用范围不变。
3. 更新两个 ED 跃迁构造和解析自由色散。根据需要选择小的有限占据参考或扩展当前 NC ED 尺寸接口；记录其基定义。
4. 从色散独立计算自由热密度、能量和动量占据。对于配对，在稳定高斯情形下使用 $`\omega_{\mathbf k}=\sqrt{(\varepsilon_{\mathbf k}-\mu)^2-\Delta^2}`$。对照对各 $`t_\nu`$ 的导数检验动能。
5. 通过直接稠密乘法验证非均匀固定 HS 场，然后运行相互作用的 BAFQMC/ED 比较。在 $`(1,1,1)`$ 处用现有论文情形恢复原始三角晶格结果。
6. 添加小型精选示例、精确命令、资源测量和模型定义到文档中。将较大的生成扫描保存在其专属输出目录中。

令 $`t_3=0`$ 会移除对角图键。要描述物理正方晶格，还需设置正交原胞基矢并更新倒格坐标和动量可观测量。这说明图的连通性与其实空间嵌入是模型的独立部分。

## 独立验证与交付物

直接从新模型定义构建小型稠密参考。它应该独立于生产辅助程序来组装格点算符、跃迁、配对和相互作用项。在共同的有限 Fock 基中比较哈密顿量能谱和热力学可观测量；还要检验厄米性、配对对称性、能量分解以及相关的热力学或规范恒等式。[test_ed_physics.py](https://github.com/YinkaiYu/BAFQMC/blob/main/tests/physics/test_ed_physics.py) 中的方法是一个起点。

对于显式非均匀场，独立地构成有序短时矩阵及其完整乘积。比较 Green 函数、包含标量在内的完整权重以及提议更新比率。测试多个场构型和稳定化间隔。将 [test_fixed_field_solver.py](https://github.com/YinkaiYu/BAFQMC/blob/main/tests/physics/test_fixed_field_solver.py) 扩展到新的几何/通道；其当前覆盖检验相互作用 HS 传播，但不测试非零 Metropolis 提议。更新公式的更改还需要针对确定性非零提议的新/旧稠密比率和更新后的 Green 矩阵。

使用 [test_gaussian_solvers.py](https://github.com/YinkaiYu/BAFQMC/blob/main/tests/physics/test_gaussian_solvers.py) 中的解析自由或 Bogoliubov 极限，然后运行一个实际包含新项的相互作用计算。不设通用统计 sigma 门槛地报告均值、SEM 和参考差异。明确测试新模型：通过旧的三角 benchmark 只证明那些情形保持一致。

在已添加的模型测试旁边运行现有检验：

```bash
make check
make physics
```

当 ED 有独立环境时使用 `PYTHON_ED=/path/to/quspin/python`。这些检验有独立的物理参考；新模型应添加同样直接的覆盖。完整的 22 点生产运行在请求的验证范围需要时，可通过 `python3 reproduce.py` 获得。

交付经过文档记录的哈密顿量和晶格、适合其参数的 HS/对称性和迹推导、可工作的求解器和参考实现、可运行的小型示例、独立测试以及测量的计算成本。将对现有输入或结果的有意更改记录为正常开发的一部分。研究人员应能够从本检出同时运行原始示例和新物理模型。
