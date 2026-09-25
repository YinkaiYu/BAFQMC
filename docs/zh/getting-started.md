# 开始使用

在 Linux 或 WSL 中运行计算工作流。以下命令均在仓库根目录下执行。
在 WSL 环境中，将代码检出到 Linux 文件系统可以获得更好的性能，
因为大量小型测量文件的写入对文件系统性能有较高要求。

## 安装环境

安装 MPI Fortran 编译器、Make 以及 BLAS/LAPACK。例如在 Ubuntu 或 Debian 上：

```bash
sudo apt-get update
sudo apt-get install -y git make gfortran openmpi-bin libopenmpi-dev libblas-dev liblapack-dev
```

使用 Conda 创建随附的 Python 环境：

```bash
conda env create -f benchmarks/paper/environment.yml
conda activate bafqmc
python3 scripts/doctor.py
```

该环境包含 NumPy、Matplotlib 以及用于 ED 的 QuSpin。
若已有 Python 3.11 环境，也可直接安装两个依赖文件：

```bash
python3 -m pip install -r benchmarks/paper/requirements.txt -r benchmarks/paper/requirements-ed.txt
```

Makefile 默认使用 `mpifort`（GNU Fortran）。已激活的 Intel MPI/MKL 环境同样支持。
构建配置与可选覆盖项详见[共享数值计算指南](https://github.com/YinkaiYu/BAFQMC/blob/main/src/common/README.md)。
若 QuSpin 使用独立的 Python 解释器，请在根复现命令和环境诊断工具中传入
`--python-ed /path/to/python`。

## 检查安装

```bash
make check
python3 reproduce.py --mode smoke
```

smoke 工作流会编译两个求解器，运行小规模蒙特卡洛计算，
计算精确参考值并检查分析流程。结果写入 `benchmarks/paper/output/smoke/`。
若需再次运行，请指定新的目标路径：

```bash
python3 reproduce.py --mode smoke --output /tmp/bafqmc-smoke-second
```

smoke 算例使用较少的采样量和 ED 基组规模，用于在正式计算前验证安装是否正确。

## 复现所有 benchmark 数据

```bash
python3 reproduce.py --plan
python3 reproduce.py --output benchmarks/paper/output/paper-run
```

第二条命令完成编译和运行 BAFQMC、计算 ED 参考值，然后从新测量数据
生成表格和图像。共涵盖 22 个参数点：8 个相对密度点、7 个配对点和 7 个总密度点。
每个生产点使用 100000 个测量 bin 以及原始输入设置。
在参考桌面 CPU 上，请预留 **12–24 小时**、**16 GiB 内存**（约 8 GiB 可供计算使用）
以及 **8 GiB 可用磁盘空间**。详见[实测资源估计](./benchmarks-paper-resources.md)。

输出目录的结构如下：

```text
runs/                                原始 BAFQMC 数据和 ED 结果
progress.json                        已完成阶段及耗时
observables.csv                      均值、SEM 和 ED 参考值
block_means.csv                      分块测量数据
records.json                         含参数和来源的结果记录
figures/benchmark_combined.pdf       相对密度和配对 benchmark 图
figures/benchmark_total_density.pdf  总密度 benchmark 图
environment.json                     环境与所选算例的元数据
```

若计算中断，使用相同的输出目录和选项继续：

```bash
python3 reproduce.py --output benchmarks/paper/output/paper-run --resume
```

已完成的阶段会被验证并复用；中断的阶段从原始输入重新开始。
根运行器自动管理临时目录；使用 `--work-dir /path/to/empty/linux-scratch`
可手动指定临时目录位置。

若只复现论文的某一部分：

```bash
python3 reproduce.py --scope combined
python3 reproduce.py --scope total_density
python3 reproduce.py --scope combined --model pairing
```

不带参数的 `python3 reproduce.py` 默认计算完整的全部算例。
`--mode plot` 可在不重新计算的情况下基于已存储的处理数据重新绘图。
所有阶段的列表和数据约定详见[复现指南](./benchmarks-paper-readme.md)。

## 开始你自己的计算

向你的 Agent 说明哈密顿量、晶格尺寸、温度、可观测量、
期望的参数扫描范围以及可用的计算资源。让它阅读
[AGENTS.md](./agents.md) 并按照
[自定义计算任务指南](./agent-workflows.md)操作。
各求解器的运行脚本接受独立的 manifest 文件，
因此研究输入可以与已发布的 benchmark 包分开开发。
