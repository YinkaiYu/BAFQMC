# 开发与验证

在修改物理核心代码之前，请先阅读 [AGENTS.md](./agents.md) 和
[算法指南](./algorithm.md)。本指南中的命令均在已安装的 Linux 环境的仓库根目录下执行。

## 定位修改位置

| 任务 | 起点文件 |
| --- | --- |
| 添加等时可观测量 | `src/*/src/obser_equal.f90`、对应的 ED 驱动程序、分析解析器 |
| 修改 HS 辅助场或局域更新提议 | `src/*/src/fields.f90`、`operator_Hubbard.f90`、`localU.f90`、`local_sweep.f90` |
| 修改稳定化传播 | `src/*/src/stabilization.f90`、`process_matrix.f90`、`multiply.f90` |
| 修改晶格或跃迁 | `src/*/src/lattice.f90`、`calc_basic.f90`、`non_interact.f90` 以及 ED 几何设置 |
| 添加或调度计算任务 | 求解器 `run_paper.py` 和 `benchmarks/campaign*.py` |
| 修改论文编排逻辑 | `benchmarks/paper/production.py`、`reproduce.py` |
| 修改统计或图像输出 | `benchmarks/paper/analysis.py`、`plot_manuscript.py` |
| 修改线性代数或随机数生成 | `src/common/` 及其数值不变量检查 |

每个求解器的 Makefile 列出了当前启用的 Fortran 源文件。特别地，
粒子数守恒求解器的 `globalK.f90` 和 `global_update.f90` 是保留的源文件，
但不在当前可执行文件中编译。在扩展代码路径之前，请先检查构建列表。

## 与修改规模相称的检查

对于文档修改，检查链接和命令，并运行 `git diff --check`。
对于 Python 或编排层修改，从快速检查开始：

```bash
make check
```

对于核心算法、估计量、数值依赖或执行流程的修改，运行与独立参考值对比的科学检查：

```bash
make physics
```

[测试指南](./testing.md)介绍了高斯 Fortran 计算、
独立有限 Fock 空间 ED、热力学恒等式检查和数值库检查。
若 QuSpin 使用独立环境，在 Make 命令中传入 `PYTHON_ED=/path/to/python`。
对于续跑/来源追踪逻辑的修改，还需运行实机流水线测试：

```bash
BAFQMC_RUN_MPI_TESTS=1 python3 -m unittest discover -s tests/reproduction -p test_production.py -v
```

若 ED 使用独立环境，请为该命令设置 `BAFQMC_PYTHON_ED=/path/to/python`。
通过公开复现入口点进行安装检查：

```bash
python3 reproduce.py --mode smoke --output /tmp/bafqmc-development-smoke
```

若需要，同时传入 `--python-ed /path/to/python`。

## 物理估计量与算法

先写下算符定义：组分求和、等时排序、归一化、连通与非连通部分，
以及适用时的配对相位。基于实际使用的 Green 函数约定推导估计量。
在 ED 和分析 schema 中添加对应项，确保从求解器到图像的每一步输出
都保有完整的物理定义。

添加一个与修改物理内容相对应的小体系直接对比。
现有的有限 Fock 空间测试在固定 Hilbert 空间中对有相互作用的 ED 进行检验；
高斯测试则将实际求解器的测量结果与未截断的解析解进行比较，
包括配对相位和可观测量归一化。

对于 HS 更新或有相互作用传播的修改，还需运行相关的
`make -C src/<solver> benchmark-dqmc` 计算，并在需要时将其扩展到修改后的参数区间。
高斯 `Nwrap` 对比覆盖二次型传播；有相互作用稳定性和严重病态的低温情形
需要在相应参数下进行专项测试。
对于随机对比，需同时记录采样误差和收敛参数。
论文复现保留所有有效结果及其不确定度。

## 保持结果可复用

在每次新的对比中记录参数、初始种子、编译器与数值库、采样量、分块设置以及 ED 截断。
跟踪测试所需的小型测试输入和处理后的参考值。
将原始链数据、日志、编译文件和本地预览存放在被 Git 忽略的输出目录中。

更新 benchmark 数据时，保持输入、参考参数、分块均值和图示统计量的一致性。
记录计算过程和修改原因，以便读者通过 Git 历史追溯。
使用 `--mode check` 重建均值和 SEM 并与 ED 参考参数进行比较。

各求解器的开发指南（含文件格式说明）：
[粒子数守恒求解器](https://github.com/YinkaiYu/BAFQMC/blob/main/docs/solvers/number_conserving/development.md) 和
[配对求解器](https://github.com/YinkaiYu/BAFQMC/blob/main/docs/solvers/pairing/development.md)。

## 双语文档

英文页面由仓库根目录 README、`docs/` 和求解器指南生成；中文页面位于
`docs/zh/`。修改任一语言时，同步相应页面的物理定义、命令和资源估计。
两套配置统一构建和部署，中文站点位于 `/BAFQMC/zh/`。
文档和图片使用相对于源文件的链接；构建脚本将其转换为站点链接，
同时保留源文件在 GitHub 上的可读性。

```bash
python3 -m pip install -r docs/requirements.txt
python3 scripts/build_docs.py
mkdir -p .build/preview
ln -sfn ../site .build/preview/BAFQMC
python3 -m http.server 8000 --directory .build/preview
```

打开 `http://localhost:8000/BAFQMC/` 或 `http://localhost:8000/BAFQMC/zh/`，
在浏览器中检查公式、语言切换、中文搜索和代码复制按钮。严格构建也会检查本地 Markdown 链接。
