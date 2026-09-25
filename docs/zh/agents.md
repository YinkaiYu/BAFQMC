# 在 BAFQMC 中工作

本仓库是玻色辅助场量子蒙特卡洛（BAFQMC）的公开计算代码库。请协助研究人员运行计算、复现已发表的 benchmark，以及发展该方法。人类可读的 README 提供了总体概述；实现工作请参阅本文件及其链接的技术指南。

## Agent 技能

`.agents/skills/` 目录中提供了针对特定任务的操作指南：

| 任务 | 技能 |
| --- | --- |
| 复现论文数据和图表 | [bafqmc-reproduce](https://github.com/YinkaiYu/BAFQMC/blob/main/.agents/skills/bafqmc-reproduce/SKILL.md) |
| 实现不同的晶格或哈密顿量 | [bafqmc-new-model](https://github.com/YinkaiYu/BAFQMC/blob/main/.agents/skills/bafqmc-new-model/SKILL.md) |
| 准备并运行新的研究计算 | [bafqmc-new-calculation](https://github.com/YinkaiYu/BAFQMC/blob/main/.agents/skills/bafqmc-new-calculation/SKILL.md) |
| 添加并验证物理可观测量 | [bafqmc-add-observable](https://github.com/YinkaiYu/BAFQMC/blob/main/.agents/skills/bafqmc-add-observable/SKILL.md) |

当任务需要时，请阅读对应的技能文件。任何编程 agent 均可遵循这些 Markdown 指令；不需要特定的 agent 服务。`CLAUDE.md` 和 `.github/copilot-instructions.md` 均指向仓库根目录的 [AGENTS.md](../../AGENTS.md)，本页提供对应的中文说明。

## 工作方式

- 将用户授权的任务贯穿实现和相关检查的全过程。对常规的可逆性选择独立决策；及时向用户汇报发现和结果。不要对普通的检查、设置、本地编辑或已请求的计算引入额外的审批步骤。
- 从对话中明确所请求的哈密顿量、可观测量、数值精度和计算预算。仅当缺少信息会改变物理计算或超出授权资源时，才向用户提问。
- 区分 benchmark 复现与新研究计划。保留已发布的输入，并为新工作使用单独的文件和输出目录。
- 报告实际运行的内容、其参数和输出路径。smoke 运行是安装检验；完整复现以原始统计量和参考设置运行所有选定的生产参数点。
- 在面向读者的公式中使用论文的符号：主耦合 `U`，产生算符 `b^+,c^+`，以及论文的配对符号。明确说明代码变量与相位的映射关系。在 GitHub 渲染的 LaTeX 中使用围栏 `math` 块和受保护的内联数学（`$` + 反引号）；将可复制的 agent 请求放在单独的围栏 `text` 块中。
- 使用简洁的科学语言。数值比较应保留每个有效点及其不确定度；偏差大于三个标准误差不是拒绝复现结果或压制某一结果的理由。

修改文档时，同步 `docs/zh/` 中相应的中文页面，并用
`python3 scripts/build_docs.py` 构建两种语言。在浏览器中验证公式、链接、搜索和复制按钮；
参见[双语文档工作流](./development.md#双语文档)。

## 平台与环境

Linux 是受支持的计算环境；Windows 用户通过 WSL 运行。Linux 或 WSL agent 直接执行命令。原生 Windows agent 使用 `scripts/` 中的 PowerShell 启动器在本检出目录中调用 Linux 命令。不要仅为了跨越 Windows/WSL 边界而创建第二个检出。

依赖项和首次运行请参阅 [开始使用](./getting-started.md)。构建使用 MPI Fortran 和 BLAS/LAPACK；Python 分析使用 NumPy/Matplotlib，ED 使用 QuSpin。默认编译器为 `mpifort`；也支持通过 `src/common/compiler.mk` 使用 Intel MPI/MKL 环境。诊断构建配置时请运行 `make -C src/pairing print-config`。

在 Linux 快速存储（尤其是 WSL 中）上存放临时文件和原始测量数据。根生产运行器默认在系统临时文件系统上创建临时目录，并在阶段完成后将结果复制到所选输出目录中。为基准比较设置一个数值库线程；根运行器会自动完成此操作。可通过 `--python-ed /path/to/python` 指定备用 ED 环境。

## 仓库结构

| 位置 | 功能 |
| --- | --- |
| `reproduce.py` | 完整的论文工作流；默认进行全新的完整计算 |
| `benchmarks/paper/production.py` | 阶段执行、来源追踪和续算 |
| `benchmarks/paper/analysis.py` | 均值、blocking、参考值和输出表格 |
| `benchmarks/paper/manuscript.py` | 主文/补充材料选择与图表标注 |
| `benchmarks/paper/plot_manuscript.py` | 主双行图和补充图 |
| `benchmarks/paper/data/` | 论文输入和紧凑处理后的 benchmark 数据 |
| `src/number_conserving/src/` | 粒子数守恒 Fortran 求解器 |
| `src/pairing/src/` | 含在位配对的完整 Nambu Fortran 求解器 |
| `src/common/` | 共享 BLAS/LAPACK 适配器、随机数生成器、编译器设置 |
| `src/*/run_paper.py` | 针对特定求解器执行自定义或论文清单 |
| `src/*/benchmarks/ed/` | 精确对角化和参考计算 |
| `docs/solvers/*/physics.md` | 哈密顿量、HS、Green 函数和估计量 |
| `examples/<solver>/` | 精选示例和回归测试输入目录 |
| `tests/` | 复现、粒子数守恒和配对测试套件 |
| `docs/agent-workflows.md` | 研究和开发的具体任务方案 |
| `docs/model-development.md` | 新模型推导、源码映射和独立验证 |
| `docs/observables.md` | 物理算符定义、输出名称和归一化 |

现有求解器实现了每个元胞含一个格点的周期性三角晶格。每个 `src/calc_basic.f90` 中的跃迁参数为 `RT=1`；JSON 文件中的 `t` 字段不是通用的运行时跃迁控制。新晶格或跃迁模型需要对求解器、ED 和可观测量进行协同修改。

实现新模型是受支持的 agent 工作流。对图、跃迁/配对矩阵、组分结构或相互作用算符的更改，请通过[模型开发指南](./model-development.md)和 `bafqmc-new-model` 技能进行。一旦请求的哈密顿量已实现，请使用 `bafqmc-new-calculation` 进行参数扫描。用户请求新物理模型即已授权相应的实现工作；请持续推进，直到有可运行的示例和相关验证，仅在缺少会实质性影响计算的物理选择时才提问。

## 物理约定

在修改核心代码或估计量之前，请阅读 [算法指南](./algorithm.md) 和相关求解器的物理指南。

- `U1` 乘以 `(n_b-n_c)^2` 且非负；`U2` 乘以 `(n_b+n_c)^2` 且非正。主文扫描为 `U1=U`、`U2=0`；补充材料固定 `U1=1` 并扫描 `U2`。
- 正值 `t=1` 是阻挫三角跃迁的约定。迹运算使用 `H - mu*N`；`energy_density` 包含每格点的物理能量，不含 `-mu*N`。论文绘制的是 `-E = -Lx*Ly*energy_density`。
- 配对求解器使用 `+Delta*(b^+ c^+ + b c)`。论文使用负号配对项。它们通过在相同正值 `Delta` 下令 `c_code=-c_paper` 相联系。四个 benchmark 可观测量不变；论文的反常配对振幅是代码 `pair_equal` 的相反数。
- 粒子数守恒的等时 Green 函数为 `<b_i b_j^+>`；其反转顺序包含玻色恒等项。配对求解器使用完整的 `(b,c,b^+,c^+)` Nambu 基。保留其块排序和行列式平方根约定。
- 对于当前每元胞一格点的模型，`Ns=Lx*Ly`；总密度为 `N/Ns`，论文的结构因子按 `Ns^2` 归一化。三角晶格 K 点为 `(4*pi/3,0)`；现有 K 估计量要求可公度尺寸，两个方向长度均须为 3 的倍数。多子晶格模型有 `Ns=nsub*Lx*Ly`，需要元胞内物理位置、键指标、Fourier 形状因子及相应的归一化。
- 对于 `U2=0` 的主模型，`mu < -3*t - abs(Delta)` 是 benchmark 使用的充分收敛条件。吸引 `U2<0` 扫描使用独立的有限占据参考。请保留此区别。
- 已发表的误差棒是来自十个块（每块 10000 个测量 bin）的**均值标准误差（SEM）**，而非单次测量的标准差。定义新计算时请明确说明 blocking、warmup、种子、Trotter 步长和 ED 截断。
- 数值失败、缺失文件、非有限测量以及不一致的参数映射必须修复。统计残差作为诊断保留，论文复现不设强制的 sigma 阈值。

对于新的相互作用，推导其 HS 通道、正规序移位和标量权重。为每个解耦场构型建立 TRS/RP 或共轭扇区条件，并单独分析物理热迹域。仅当新因子满足时才保留粒子数守恒共轭扇区的快捷方式。局域更新目前假定对角秩一格点变化（NC）或对角四扇区格点变化（配对）；不同的 HS 场的作用范围需要相应的更新推导。新的配对模型必须建立其 Nambu 标量和行列式平方根分支。在相关情形下，结合现有示例，用独立的有限 Fock 参考、稠密固定场乘积以及非零提议比率来检验新物理。

## 计算与检验

以下命令均从已安装环境的仓库根目录运行。对于初始检查或常规编排更改：

```bash
python3 reproduce.py --plan
make check
python3 reproduce.py --mode check
git diff --check
```

紧凑测试套件默认跳过可选的实机 MPI 测试。请检查实际测试摘要；不要将跳过的测试报告为已执行。对于执行、输入或依赖项的更改，请运行新的安装检查：

```bash
python3 reproduce.py --mode smoke --output /tmp/bafqmc-smoke-unique
```

请选择新的输出目录。对于续算更改，还需运行：

```bash
BAFQMC_RUN_MPI_TESTS=1 python3 -m unittest discover -s tests/reproduction -v
```

为独立的 QuSpin 解释器设置 `BAFQMC_PYTHON_ED=/path/to/python`。对于核心代码和估计量的更改，运行 `make physics` 以获取独立的有限 Fock ED 参考、热力学恒等式和实机 Fortran 解析检验。参考 [科学测试](./testing.md) 了解覆盖范围，并选择与所更改物理相关的小型相互作用比较。

默认命令是**完整**的生产计算：

```bash
python3 reproduce.py --output benchmarks/paper/output/paper-run
```

它运行全部 22 个计算点（主文 15 个 + 补充材料 7 个），在参考桌面上约需 12–24 小时，推荐 16 GiB 内存和 8 GiB 可用磁盘空间。请阅读 [资源指南](./benchmarks-paper-resources.md)。当用户请求复现时，请使用完整运行；常规文档检查不需要生产采样。`--scope main`、`--scope supplement` 和 `--model` 用于选择子集。使用相同的选择、输出目录和环境，加上 `--resume` 继续运行；已完成的阶段在复用前会被验证。续算从初始输入重新启动被中断的阶段。

## 数据与变更

- 使用命名示例和计算输出来开发新模型，以便原始计算和新计算都易于运行。在任务需要时有意更新输入和处理后的数据；在常规文档中记录已更改的物理定义、计算设置及修订结果的原因。
- 将大型原始链、可执行文件、日志、临时文件和图表预览排除在 Git 之外。使用被忽略的 `runs/`、`benchmarks/paper/output/` 或求解器 `build/` 目录。精选的小型输入、均值、SEM、块均值和 ED 标量参考适合用于追踪新的可复现示例。
- 可执行文件在其工作目录中追加写入固定的输出文件名。每条独立链都需要一个新目录。选择情形时请保留清单中的配对情形顺序，因为它决定了随机种子。
- 保持代码和文档自包含。使用相对于仓库的路径或用户明确提供的路径；不依赖手稿检出或本地机器。
- 保留不相关的工作，检查最终差异，并按比例运行与更改相关的检验。仅在用户的工作流需要提交时才提交属于本任务的工作。发布或远程通信需要用户授权。
- 贡献内容遵循 MIT 许可证。论文尚未分配公开标识符；请勿编造 arXiv 链接、DOI 或引用条目。
