# Benchmark 复现

默认的 `python3 reproduce.py` 会重新运行全部 22 个 BAFQMC/参考计算，
处理新的测量数据，并绘制随附研究的主图和补充材料 benchmark 图。
请参阅 [RESOURCES.md](./benchmarks-paper-resources.md)，了解 12–24 小时的生产预算、
内存、临时存储以及续跑命令。可选的 `--mode plot` 可用已存储的小型数据集重新绘图。
两条路径均绘制合并后的主图 benchmark 和补充材料中单独的吸引密度 benchmark。

## 数据与参考约定

| 图/行 | 变化参数 | 固定参数 | 参考 |
|---|---|---|---|
| 主图 Fig. 2(a–d) | U = U2 = 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2 | U1 = 0, mu = -3.5, beta = 4 | 粒子壳层 ED；U = 0 时的解析自由玻色子参考值 |
| 主图 Fig. 2(e–h) | Delta = 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3 | U1 = 0, U = U2 = 1, mu = -5, beta = 4 | 稠密 ED，nmax = 3，ncut = 4 |
| 补充材料 Fig. S1(a–d) | U1 = -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0 | U2 = 1, mu = -7, beta = 1 | U1 < 0 时的低密度有限占据窗口参考值 |

存储的输入和数据键保留 `U1` 和 `U2`。只有主图文字使用 `U` 代指 `U2`；
补充材料模型保留两个相互作用系数。处理后的表格和图像元数据明确记录了这一命名约定。

`--scope main` 选取主图中的 15 个点，`--scope supplement` 选取 7 个吸引密度点。
默认的 `--scope all` 包含全部。范围选择和求解器选择相互独立：

| 选择 | 点数 | 扫描 |
|---|---|---|
| `--scope main` | 15 | U 和 Delta |
| `--scope main --model number_conserving` | 8 | U |
| `--scope main --model pairing` | 7 | Delta |
| `--scope supplement` | 7 | U1 |
| `--model number_conserving` | 15 | U 和 U1 |

主图模型在每个绘图点均满足 `mu < -3*t - abs(Delta)`。
论文补充材料"辅助场域内的收敛性"子节已证明，其二次型传播和幺正相对密度 HS 因子
在整个辅助场域内给出有限的迹。该结论适用于 `U1=0` 的排斥模型；
补充材料扫描使用独立的双通道模型及其指定的有限占据参考值。

配对求解器写出实数配对项时系数为正，而手稿写作 `-Delta (b^+ c^+ + b c)`。
相位约定 `c_code = -c_paper` 在同一正扫描值下将两者联系起来。
密度、总物理能量和两个绘制的结构因子在该变换下不变。
求解器的反常配对振幅 `pair_equal` 以手稿算符表示时符号相反；
详见[配对约定](https://github.com/YinkaiYu/BAFQMC/blob/main/docs/solvers/pairing/physics.md)。

所有点的参数：t = 1，周期边界 3×3 三角晶格，Delta tau = 0.01，
100000 个测量 bin，1 个 MPI 进程。误差棒为 10 个连续分块（每块 10000 个 bin）的标准误差。
预热和稳定化设置、提议宽度以及精确初始种子均在各 `inputs/<case>/` 目录中。
配对求解器的 `seeds.txt` 包含每进程种子列表；论文使用其第一个条目，
而不是某个随机数流的多字状态。

四个测量文件为 `density_total`、`energy_density`、`sf_K` 和 `dw_K`。
结构因子文件存储实部和虚部。图像使用其实部，
并以 -E = -9e 的形式展示物理总能量，能量标准误差乘以 9。
能量中不含化学势项。复现过程中不设统计接受准则来移除数据点。

负 U1 补充参考值使用原始低密度占据窗口：6 个已完成的粒子壳层，
配置最大值为 8。这是由 ED 记录中 `low_density_cutoff_accepted` 编码的有限窗口对比。
主图相互作用扫描中的 U = 0 参考值由解析玻色分布直接计算；
其原始截断 ED JSON 仅为来源记录保留，不用于该红色数据点。
所有其他参考值从存档的 ED 结果中读取。

## 目录结构

```text
data/index.json                  所选算例、参数、绘图值和来源信息
data/observables.csv             紧凑型处理后均值、SEM 和 ED 值
data/block_means.csv             每个算例和可观测量的十个分块均值
data/<model>/manifest.json       可执行的论文计算方案
data/<model>/inputs/<case>/      精确的初始输入文件和 ED 参数
data/<model>/ed/<case>.json      原始 ED 结果
data/<model>/raw/<case>.tar.gz   可选的本地原始存档（Git 忽略）
analysis.py                      原始数据统计、参考值、表格、绘图调度
plot_manuscript.py               主 U/Delta 网格和补充 U1 扫描
plot_number.py                   原始 2×4 粒子数守恒布局
plot_pairing.py                  原始 1×4 配对布局
plot_style.py                    共享的原始图像样式
```

已跟踪的数据包约 250 kB。大型原始测量文件由 BAFQMC 生成，不进入 Git。
小型 ED JSON 文件包含参考期望值和截断元数据。存储的分块均值允许
在不传输完整蒙特卡洛链的情况下独立重建绘图均值和标准误差。
`--mode check` 无需生成图像即可检查有限记录、ED 参数以及均值和 SEM 的重建。
`--mode raw` 在可选的原始存档存在时从中重建统计量。

当前输出为 `figures/benchmark_combined.pdf`（两行）和
`figures/benchmark_attractive.pdf`（一行），旁边附有 PNG 预览。
只绘制所选范围和求解器覆盖的图像和行。
单独的主图行写作 `figures/benchmark_combined_number_conserving.pdf` 或
`figures/benchmark_combined_pairing.pdf`，分别保留面板标签（a–d）或（e–h）。
原始布局模块作为可复用的绘图辅助工具保留。

生成的表格包含 `benchmark`、`section`、`figure`、`figure_row` 和 `panels`
用于定位当前手稿中的每个点。`scan_parameter` 是绘图符号，
`solver_parameter` 是其输入键，`U` 给出主图相互作用系数。
JSON 记录在 `manuscript` 字段下汇总这些内容。
现有算例 ID 和存档的 `row` 字段保留其历史含义；
当前图像排列请使用 `figure_row`。

## 重新运行计算

`--mode full`（默认）编译每个求解器，执行论文计算方案，运行 ED，
对新测量结果进行分块，并从新结果生成两张图。原始数据保留在 `data/` 中不变。
输出写入 `--output` 指定的目录，并包含与 ED 的数值差异。

新生成的蒙特卡洛轨迹可能因编译器和数值库而异；
已发布的数据仍是重现原始曲线的精确来源。
PDF 元数据和字体库版本也可能改变 PDF 字节哈希，
即使绘图数据和布局完全相同。

独立阶段：`--mode dqmc` 和 `--mode ed` 写入 `output/runs/`。
两个阶段完成后，`--mode analyze` 对输出分块并绘制图像。
每个阶段以及 `--resume` 请使用相同的 `--scope` 和 `--model`。
各求解器的 `run_paper.py --help` 还提供算例选择或方案配置以及分析入口点。
创建新模型方案时不要覆盖存档的输入文件；请将其复制到自己的配置目录中。

## 来源记录

粒子数守恒实现从 `code_bosonDQMC` 的提交
`68b82365ab817fd4a96e357358c31181f7acb3f3` 导入；
配对实现从 `code_bosonDQMC_paring` 的提交
`e4c1d04e6b756eb6bd6d0cf5b9a778c6b54a920c` 导入。
两者的 Fortran 模型核心代码未经改动直接导入。
本发行版提供可移植的构建文件、本地计算方案运行脚本和标准数值接口；
详见[数值支持指南](https://github.com/YinkaiYu/BAFQMC/blob/main/src/common/README.md)。

这些标识符描述来源历史；复现所需的所有源文件、计算方案输入、
处理后数据和参考实现均包含在本仓库中。
`data/index.json` 记录了所选算例、物理参数和数据来源。
生产计算直接从包含的输入文件生成新的原始测量数据。
