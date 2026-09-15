# 智能体研究工作流程

首先阅读 [AGENTS.md](./agents.md)。此处的命令在 Linux 或 WSL 环境下从仓库根目录运行。本仓库提供可用的三角晶格求解器和可编辑的计算方案输入，以及用于实现不同物理模型的源码级工作流程。

`.agents/skills/` 中的可移植任务技能负责路由这些工作流程：
[复现 benchmark](https://github.com/YinkaiYu/BAFQMC/blob/main/.agents/skills/bafqmc-reproduce/SKILL.md)、
[实现新模型](https://github.com/YinkaiYu/BAFQMC/blob/main/.agents/skills/bafqmc-new-model/SKILL.md)、
[运行新计算](https://github.com/YinkaiYu/BAFQMC/blob/main/.agents/skills/bafqmc-new-calculation/SKILL.md)，以及
[添加可观测量](https://github.com/YinkaiYu/BAFQMC/blob/main/.agents/skills/bafqmc-add-observable/SKILL.md)。

## 设置并建立可运行的基准环境

1. 按照[环境配置](./getting-started.md)操作，然后运行
   `python3 scripts/doctor.py` 和 `make check`。
2. 运行 `python3 reproduce.py --mode smoke`，检验小规模 BAFQMC、ED 与分析流程。
3. 检查运行日志和结果表格，将执行的计算案例和输出目录（包括跳过了哪些可选检查）报告给用户。

若使用独立的分析和 ED 环境：

```bash
python3 scripts/doctor.py --python-ed /path/to/quspin/python
python3 reproduce.py --mode smoke --python-ed /path/to/quspin/python
```

## 实现不同的晶格或哈密顿量

当所需的晶格图、跃迁、配对、组分结构或相互作用需要修改源码时，使用 [bafqmc-new-model](../../.agents/skills/bafqmc-new-model/SKILL.md)。[模型开发指南](./model-development.md)将推导过程映射到实际的 Fortran 和 ED 符号，并给出各向异性三角晶格的实现蓝图。以下请求可直接复制给智能体：

```text
阅读 AGENTS.md 并使用 bafqmc-new-model。实现具有独立实数跃迁振幅 t1、t2、t3
及现有两个在位相互作用通道的各向异性三角晶格模型，保留各向同性模型作为示例。
推导对称性与迹收敛条件，更新求解器输入、ED 和能量估计量；用独立自由能谱、
固定场矩阵乘积和小规模相互作用计算验证，提供运行命令和实测资源消耗。
```

```text
使用 bafqmc-new-model 在 kagome 晶格上实现双组分相对密度相互作用
U (n_b - n_c)^2。两个组分采用相同的正实数最近邻跃迁，不含配对项。
明确定义三个子晶格和物理格点数，建立 HS 解耦后的权重非负条件，
实现配套的 BAFQMC 与有限 Fock 空间 ED，交付经过验证的小型示例，
并写清观测量定义和计算资源需求。
```

用户可以将上述模型定义替换为自己的定义。请说明物理哈密顿量、几何结构和边界条件、目标可观测量以及可用计算资源，智能体将提供实现方案并完成代码编写和验证。目前暂不支持通用的 `lattice` JSON 标志。添加模型需要同时实现新的输入格式及读取这些输入的数值核心。

对于多子格情形，需将物理格点与元胞及 Nambu 扇区分开处理。对于新的 HS 相互作用，既要推导场因子，也要推导相应的局域更新及其作用范围。对于配对情形，需保持标量权重、分支和物理缩并与所选 Nambu 约定的一致性。已发表的论文测试应保持可运行，新的独立测试必须覆盖新模型。实现完成后，常规扫描可使用以下计算方案。

## 复现完整论文结果

运行 `python3 reproduce.py --plan`，阅读 [RESOURCES.md](./benchmarks-paper-resources.md)，并告知用户计算预算。当用户要求完整复现时，执行以下计算：

```bash
python3 reproduce.py --output benchmarks/paper/output/paper-run
```

默认包含全部 22 个生产参数点、新生成的参考值、统计处理和图形。保留原始采样和 ED 设置。若执行中断，可用以下命令继续：

```bash
python3 reproduce.py --output benchmarks/paper/output/paper-run --resume
```

使用相同的范围、模型选择、解释器和数值库线程设置。已完成的阶段在重用前会经过验证；部分完成的阶段从初始输入重新启动。在报告完成前，请检查 `progress.json`、`observables.csv` 和两个生成的图形。

## 准备一个粒子数守恒研究计算点

使用独立的 manifest 和输入文件。以下脚本将小型自由玻色子示例复制到一个被忽略的研究工作目录中：

```bash
python3 - <<'PY'
import json
from pathlib import Path
import shutil

root = Path.cwd()
source = root / 'examples/number_conserving/examples/triangle_3x3_free'
work = root / 'runs/custom-number'
inputs = work / 'inputs/my-point'
inputs.mkdir(parents=True)
for name in ('paramC_sets.txt', 'confin.txt', 'seeds.txt', 'params.json'):
    shutil.copy2(source / name, inputs / name)
parameters = json.loads((inputs / 'params.json').read_text())['parameters']
manifest = {
    'cases': [{
        'id': 'my-point',
        'input_dir': 'inputs/my-point',
        'ed_params': 'inputs/my-point/params.json',
        'parameters': parameters,
        'analysis': {'block_size': 2, 'skip_samples': 0},
    }]
}
(work / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
PY
```

这是一个包含八个测量 bin 的安装检查示例。若要修改物理设置，需同步编辑 `paramC_sets.txt`、`params.json` 中的 `parameters` 对象以及 manifest 中该计算案例的 `parameters`。第一个数字行为 `U1 U2 mu`，第二行为 `Lx Ly Ltrot beta`。第四行包含 `Nwrap Nbin Nsweep shiftLoc`；`dtau=beta/Ltrot`。在超出安装示例范围时，应复制适当的正式计算预热设置。ED 的粒子壳层和收敛策略在 `params.json` 中。选择合适的块大小，使块数足以根据所需采样数估计 SEM。

检查命令后，对同一输出目录依次执行以下三个阶段：

```bash
python3 src/number_conserving/run_paper.py --manifest runs/custom-number/manifest.json --output runs/custom-number/results --mode dqmc --dry-run
python3 src/number_conserving/run_paper.py --manifest runs/custom-number/manifest.json --output runs/custom-number/results --mode dqmc
python3 src/number_conserving/run_paper.py --manifest runs/custom-number/manifest.json --output runs/custom-number/results --mode ed
python3 src/number_conserving/run_paper.py --manifest runs/custom-number/manifest.json --output runs/custom-number/results --mode analyze
```

在独立环境中运行 ED 时，添加 `--python /path/to/quspin/python`。一个 manifest 可包含多个各自拥有独立输入的计算案例。`--case ID` 可重复使用以选择子集。输出包括 `my-point/dqmc/`、`my-point/ed/results.json` 和 `analysis.json`。

对于不含 ED 的蒙特卡洛计算，也可使用直接输入目录运行器；参见[粒子数守恒求解器指南](./solver-number-conserving.md)。

## 准备配对研究扫描

配对 manifest 是自包含的：其参数会生成求解器和 ED 的输入文件。复制小型流水线示例：

```bash
mkdir -p runs/custom-pairing
cp src/pairing/benchmarks/campaigns/triangle_pairing_pipeline_smoke.json runs/custom-pairing/manifest.json
```

编辑 `parameters` 中的 `Lx`、`Ly`、`U1`、`U2`、`mu` 和 `beta`；编辑 `deltas` 进行扫描。当前跃迁参数保持 `t=1`。`dqmc_defaults` 包含 `dtau`、测量块数、预热步数、提议宽度和分块方式；`ed_defaults` 包含 `nmax` 和 `ncut`。根据所需分析一致性设置 `analysis_defaults.block_size`。复制的示例从八个测量 bin和较小的 ED 基开始。已发表的[配对 manifest](https://github.com/YinkaiYu/BAFQMC/blob/main/benchmarks/paper/data/pairing/manifest.json) 提供完整采样设置作为参考。

```bash
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode init --seed 24681357 --dry-run
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode init --seed 24681357
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode dqmc --seed 24681357
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode ed --seed 24681357
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode analyze --seed 24681357
```

若在独立 ED 解释器中运行，添加 `--python /path/to/quspin/python`。生成的计算案例名称在试运行中显示。使用可重复的 `--case` 参数选择案例，同时保留完整 manifest 及其顺序：顺序决定每个案例的初始随机种子。输出包括 `inputs/<case>/`、`ed_results/<case>.json` 和 `summary/comparison_observables.csv`。

在扩展所需研究扫描之前，先运行一个规模适中的案例，测量时间和内存消耗。为独立链选择新的输出目录。底层可执行文件将测量结果追加到固定文件名中。

## 添加可观测量

1. 写出其算符定义和单位。说明等时顺序、组分求和、归一化、动量以及是否需要做关联减法。
2. 阅读求解器的 `docs/solvers/<solver>/physics.md`，推导 Wick 估计量。保留 `<b_i b_j^+>` 中的玻色子恒等贡献，以及配对情形中完整的 Nambu 块约定。
3. 扩展 `src/obser_equal.f90` 及其输出、对应的 ED 可观测量，以及分析解析器/模式。将定义放在实现旁边。
4. 用解析极限、ED 恒等式或独立 ED 表示验证一个小系统案例，然后运行适当的实机 BAFQMC 对比。
5. 提供可复现的输入，并报告均值、SEM 和参考值。
   在有用的情况下追踪紧凑的测试样例，将大型原始输出保存在本地。

使用[开发指南](./development.md)获取精确的检查命令。对于哈密顿量本身的修改，继续使用[新模型工作流程](./model-development.md)，包括 HS 构造、对称条件、几何结构、ED 构造和模型特定测试。
