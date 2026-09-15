# 粒子数守恒 BAFQMC 与精确对角化

`src/number_conserving/` 包含用于论文中双组分三角晶格 Bose-Hubbard 哈密顿量的 Fortran 有限温度 BAFQMC 求解器、QuSpin ED 实现、回归测试用例以及计算分析工具。哈密顿量、Hubbard–Stratonovich 场、Green 函数、估计量及截断条件详见 [物理指南](https://github.com/YinkaiYu/BAFQMC/blob/main/docs/solvers/number_conserving/physics.md)。

该求解器产生正文 benchmark 中 `U1=0`、`U2=U` 的相互作用扫描结果，以及补充材料中在 `U2=1` 下改变 `U1` 的吸引密度扫描结果，两者均有 `Delta=0`。输入名称 `U1` 和 `U2` 在整个实现中分别标识总密度和相对密度通道。

## 编译与运行

在 Linux/WSL 上安装 MPI Fortran 编译器和 BLAS/LAPACK；默认编译支持使用 `mpifort` 的 GNU Fortran。Python 分析需要 NumPy，绘图需要 Matplotlib。ED 还需要 QuSpin、SciPy 和 Numba。依赖安装详见[安装指南](./getting-started.md)。

```bash
make -C src/number_conserving build
make -C src/number_conserving run-example
make -C src/number_conserving benchmark-fast
make -C src/number_conserving check-fixtures
```

`run-example` 将八个 bin 的 3x3 自由玻色子示例复制到 `src/number_conserving/build/example` 后运行一个 MPI 进程。若输出目录非空则拒绝运行。使用 `RUN_DIR=/path/to/new/run` 指定其他输出位置。若要开发新的模型参数点，将 `examples/number_conserving/examples/triangle_3x3_free` 中的三个输入文件复制到新运行目录，编辑 `paramC_sets.txt`，然后运行：

```bash
bash src/number_conserving/scripts/run_local.sh /path/to/new/run 1
```

可执行文件从工作目录读取 `paramC_sets.txt`、`confin.txt` 和 `seeds.txt`，并将测量结果追加写入该目录，同时写入重启状态。开始新链时请勿在旧输出目录中重新运行。

物理计算核心源代码保留自 `code_bosonDQMC` 提交 `68b82365ab817fd4a96e357358c31181f7acb3f3`。[共享数值适配器](https://github.com/YinkaiYu/BAFQMC/blob/main/src/common/README.md) 提供可移植编译所使用的 BLAS/LAPACK 接口。处理后的 benchmark 数据及其初始输入位于 `benchmarks/paper/data/` 下。

## 论文生产入口

根目录的复现命令负责准备论文 manifest。该模式专用入口也可独立运行任意选定的用例：

```bash
python src/number_conserving/run_paper.py \
  --manifest /path/to/manifest.json --output /path/to/new/results --mode dqmc
python src/number_conserving/run_paper.py \
  --manifest /path/to/manifest.json --output /path/to/new/results --mode ed
python src/number_conserving/run_paper.py \
  --manifest /path/to/manifest.json --output /path/to/new/results --mode analyze
```

在仓库根目录执行 `python3 reproduce.py --scope main --model number_conserving` 可复现正文的八个 U 数据点，执行 `python3 reproduce.py --scope supplement` 可复现补充材料的七个 U1 数据点。其物理参数和采样参数列于[复现指南](./benchmarks-paper-readme.md)。

使用支持 QuSpin 的解释器运行 ED，或传入 `--python /path/to/python`。使用 `--case ID`（可重复）选择用例，使用 `--dry-run` 预览命令。Manifest 中的路径均相对于 manifest 文件：

```json
{
  "analysis": {"block_size": 1000, "skip_samples": 0},
  "cases": [{
    "id": "example",
    "input_dir": "inputs/example",
    "ed_params": "inputs/example/params.json",
    "parameters": {
      "Lx": 3, "Ly": 3, "t": 1.0,
      "U1": 0.0, "U2": 1.5, "beta": 4.0, "mu": -3.5
    }
  }]
}
```

输出位于 `OUTPUT/ID/dqmc/`、`OUTPUT/ID/ed/results.json` 和 `OUTPUT/analysis.json`。原始 ED 参数和收敛策略保持不变读取。对于稳定的自由玻色子参数点，`ed` 阶段会生成新的精确解析参考值，与论文处方一致。分析阶段导出所有有效样本的均值、标准误差和与 ED 的差值。缺失或格式错误的输出视为错误。

完整的 3x2 实机回归套件为 `make -C src/number_conserving benchmark-dqmc`；它使用原始的含 100000 个 bin 的相互作用用例，耗时约 10–15 分钟。该套件独立于较长的 3x3 论文生产计算。

## 源代码结构

以下路径相对于 `src/number_conserving/`。

- `src/`：模型、HS 场、局域更新、稳定传播、可观测量。
- `benchmarks/ed/EDtriangle_quspin_3x3.py`：论文 ED，使用粒子数和平移对称块、壳层收敛诊断及增量结果。
- `benchmarks/ed/EDtriangle_symm_NEblock.py`：传统 3x2 ED 参考驱动程序。
- `benchmarks/campaign_analysis.py`：原始解析器、分块统计、归一化、精确自由参考值和可靠性诊断。
- `benchmarks/campaign_3x3.py`：通用计算输入/报告辅助工具；其默认计算范围比所选论文图更宽。
- `../../examples/number_conserving/`：精简示例、smoke 测试和实机回归输入。
- `../../tests/number_conserving/`：Python 单元测试和可选运行时测试。
- `../../docs/solvers/number_conserving/`：物理和开发指南。
- `benchmarks/fixtures/`、`references/`、`dqmc_references/`：回归数据。
- `../common/`：共享可移植数值适配器。

绘图工作流及归档的论文输入/数据记录于仓库的复现入口。

若要快速执行完整的流程检查，将 `--manifest src/number_conserving/benchmarks/campaigns/pipeline_smoke.json` 传入上述三个模式。该 manifest 运行八个自由玻色子 bin 并计算精确自由参考值。可选的 QuSpin 运行时测试单独执行一个微型相互作用 ED 计算。历史上的含相互作用 3x2 示例单独保留，不作为默认入门运行。

验证：五用例实机回归套件在可移植 Intel/MKL 编译版本上通过（所有四个相互作用用例均保留 100000 个 bin）。Python 测试套件、可选的两个 MPI 运行时测试、微型 QuSpin 运行时测试以及内存保护检查点测试均已通过。

ED 驱动程序和 `run_paper.py` 接受 `--dense-memory-cap-gib` 参数（默认为 12）。在进行稠密对角化之前，会使用保守的四复数矩阵估算进行检查；若超出预算，则在保留最后一个已完成粒子壳层的同时停止计算。此保护机制不改变物理截断条件。归档的最大相互作用论文分块维度为 9075（第七壳层），稠密工作空间估算为 4.91 GiB；第八壳层维度将达到 27225（44.18 GiB）。串行完整复现建议使用至少 16 GiB 物理内存、约 8 GiB 可用内存的机器。
