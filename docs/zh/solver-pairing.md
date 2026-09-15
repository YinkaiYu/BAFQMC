# 配对 BAFQMC 与精确对角化

本求解器对含格点配对的双组分三角晶格 Bose–Hubbard 模型进行采样。其完整四扇区 Nambu Green 矩阵和行列式平方根实现了论文中所用的配对构造。Fortran 物理计算核心保留自 `code_bosonDQMC_paring` 提交 `e4c1d04e6b756eb6bd6d0cf5b9a778c6b54a920c` 的源代码快照。

以下命令在 Linux 或 WSL 中从 `src/pairing/` 目录运行。环境配置详见仓库[安装指南](./getting-started.md)，完整论文复现命令见根目录 README。

本求解器产生联合正文 benchmark 的面板（e–h）：在 `U=1`、`mu=-5`、`beta=4` 下的七个 Delta 数据点。在求解器输入中设置 `U1=0` 和 `U2=U`；这些名称保留了总密度和相对密度 HS 通道。求解器的正实配对系数与论文中的负配对项的关系为 `c_code=-c_paper`。所有四个绘图可观测量不受此相位约定影响；反常振幅 `pair_equal` 在论文算符中符号相反。详见 [哈密顿量和相位约定](https://github.com/YinkaiYu/BAFQMC/blob/main/docs/solvers/pairing/physics.md)。

## 编译与首次运行

```bash
make build
make run-example
```

编译使用共享的 `../common/` LAPACK 适配器和随机数生成器，需要 MPI Fortran 编译器和 BLAS/LAPACK。激活 Intel MPI/MKL 环境后同样支持。如有需要，可通过 Make 覆盖 `FC`、`FFLAGS`、`LDFLAGS` 或 `LDLIBS`。

`make run-example` 将已提交的输入复制到 `build/example/`，然后运行一个 MPI 进程。后续运行请选择另一个空目录：

```bash
make run-example RUN_DIR=build/my-example MPI_NP=1
```

修改复制后的 `paramC_sets.txt` 以探索其他模型或蒙特卡洛参数；其七行数值记录于[开发指南](https://github.com/YinkaiYu/BAFQMC/blob/main/docs/solvers/pairing/development.md)。跃迁系数在 `src/calc_basic.f90` 中为 `RT=1`。运行输出使用工作目录中的固定文件名。可执行文件追加写入可观测量，因此每条独立链需要一个新目录。`seeds.txt` 每个 MPI 进程包含一个标量种子。

## 论文计算与分析

`run_paper.py` 接受所附的论文 manifest 或使用相同模式的自定义 manifest。各阶段共用同一输出目录：

```bash
python run_paper.py --manifest ../../benchmarks/paper/data/pairing/manifest.json --output build/paper --mode init
python run_paper.py --manifest ../../benchmarks/paper/data/pairing/manifest.json --output build/paper --mode dqmc
python run_paper.py --manifest ../../benchmarks/paper/data/pairing/manifest.json --output build/paper --mode ed
python run_paper.py --manifest ../../benchmarks/paper/data/pairing/manifest.json --output build/paper --mode analyze
```

生产计算包含七个配对值，每个值 100000 个测量 bin。这些是生产级计算；快速安装检查请使用仓库 smoke 命令。在仓库根目录执行 `python3 reproduce.py --scope main --model pairing` 可运行此扫描并生成带有论文面板标签（e–h）的 `benchmark_combined_pairing.pdf`。默认种子策略为归档计算策略：从 `base_seed + case_index` 初始化每个用例并生成其每进程种子列表。使用 `--seed` 启动独立链，使用 `--np` 设置 MPI 进程数。结果随编译器、数值库和进程数而变化。使用 `--dry-run` 在不写入文件的情况下检查参数和每用例种子。重复使用 `--case CASE_NAME` 执行选定用例；初始化始终使用完整 manifest，为单独调度的用例保留原始种子位置。在 WSL 中，请选择 Linux 文件系统上的输出目录而非 `/mnt/c`：求解器会追加写入大量小的可观测量记录。

在 AMD Ryzen 5 9600X（WSL2 环境、Intel Fortran/MKL、一个 MPI 进程、一个数值库线程）上的标定结果：在论文物理参数（含原始 500 次热化迭代）下，1000 个 bin 耗时 46.38 秒，峰值 RSS 为 65 MiB。线性外推得全部七个 100000-bin 运行约需九小时。原始归档运行记录的累计用例耗时为 12.2 小时。在类似机器上，配对 BAFQMC 计算约需 9–13 小时，具体取决于系统负载和存储速度。其完整可观测量输出约占 683 MiB。

各阶段写入 `inputs/<case>/`、`ed_results/<case>.json` 和 `summary/comparison_observables.csv`。JSON 比较还保留了分块误差、虚部诊断和对截断敏感的可观测量。其历史文件名 `comparison_dqmc_ed_nmax3_ncut4.json` 与论文绘图输入一致；自定义运行的截断值请查阅嵌入的 manifest。分析阶段导出所有有效均值、标准误差和参考差值。`--require-agreement` 可选地启用继承的三标准误差检验。

能量文件存储每格点的物理能量，不含化学势项。论文绘制的是带负号的总物理能量：`-E = -Lx * Ly * energy_density`。

## ED 与验证

激活仓库的 QuSpin 环境，或向 `run_paper.py` 传入 `--python /path/to/python`。独立运行 ED：

```bash
python benchmarks/ed/ed_pairing_triangle_general.py --params benchmarks/ed/params_triangle_pairing_smoke.json --output build/ed-smoke.json
make benchmark-ed
make check-fixtures
make benchmark-fast
```

`make benchmark-ed` 检验两种独立的微型基底 ED 表示及测量可观测量的有限差分恒等式。`make benchmark-fast` 运行一个小规模的实机 BAFQMC 比较。`make benchmark-dqmc` 还将零配对极限与内置的粒子数守恒参考结果进行比较。

配对 ED 在规定的占据截断范围内执行完整热迹计算。论文使用 `nmax=3,ncut=4`（7297 个基矢），稠密对角化峰值估算约 3.2 GiB。用例串行运行；默认内存上限为 12 GiB。仅在与可用内存匹配时增大 `--dense-memory-cap-gib`。更大的截断值增长迅速，属于独立的科学计算。在同一台机器和单线程下，`Delta=0.2` 的一个实际论文规模 ED 用例耗时 40.40 秒，峰值 RSS 为 3.01 GiB；其四个论文可观测量与归档 ED 值的绝对偏差在 `5.3e-18` 以内。七个 ED 用例在上述条件下约需五分钟，此时间在 BAFQMC 时间之外。上述计时于 2026-09-15 测定。

哈密顿量、Nambu 约定、HS 标量因子、Green 函数块和 Wick 估计量详见[物理指南](https://github.com/YinkaiYu/BAFQMC/blob/main/docs/solvers/pairing/physics.md)；与粒子数守恒求解器的可观测量对应关系由[可观测量约定](https://github.com/YinkaiYu/BAFQMC/blob/main/docs/solvers/pairing/pairing_observable_contract.md)定义。
