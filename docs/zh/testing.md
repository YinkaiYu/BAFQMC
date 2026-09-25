# 科学测试

测试程序将数值结果与独立计算的物理参考值进行对比。
[完成环境配置](./getting-started.md)后，在仓库根目录运行：

```bash
make check
make physics
```

`make check` 检查计算输入、数据来源、可观测量归一化、分块和分析流程。
`make physics` 编译数值代码，并将物理结果与独立参考值进行对比。
使用独立 Python 环境时：

```bash
make check PYTHON=python3 PYTHON_ED=/path/to/quspin/python
make physics PYTHON=python3 PYTHON_ED=/path/to/quspin/python
```

## 物理检查的内容

| 检查项 | 独立参考 | 检验的物理量 |
| --- | --- | --- |
| [有限 Fock ED](https://github.com/YinkaiYu/BAFQMC/blob/main/tests/physics/test_ed_physics.py) | NumPy 升降算符矩阵和在 81 态两格点 Hilbert 空间中的直接密度矩阵迹 | 两种配对 ED 实现：谱、自由能、密度、动能和相互作用能、配对振幅、配对结构因子以及粒子数矩 |
| [高斯 Fortran 求解器](https://github.com/YinkaiYu/BAFQMC/blob/main/tests/physics/test_gaussian_solvers.py) | 未截断二次型模型的动量空间 Bogoliubov 解析解 | 两种 BAFQMC 求解器：论文的全部四个可观测量、反常振幅、零配对极限、配对相位变换以及稳定化区间 |
| [固定有相互作用 HS 场](https://github.com/YinkaiYu/BAFQMC/blob/main/tests/physics/test_fixed_field_solver.py) | 对给定非均匀场直接用 NumPy 计算时间片乘积、逆矩阵和行列式 | 粒子数守恒的权重、Green 函数诊断量、格点密度、动能与相互作用能，以及两个虚时扫描方向 |
| [ED 可观测量恒等式](https://github.com/YinkaiYu/BAFQMC/blob/main/src/pairing/benchmarks/ed/test_observable_identities.py) | 热力学导数和算符恒等式 | 从自由能配对导数得到的配对振幅、从跃迁导数得到的动能、能量分解和玻色粒子数归一化 |
| [数值接口](https://github.com/YinkaiYu/BAFQMC/blob/main/src/common/README.md) | 矩阵重建和已知随机数递推关系 | BLAS/LAPACK 运算以及两个求解器共用的随机数生成器 |

有限 Fock 参考程序自行构建哈密顿量、算符和热力学迹，并在两种生产 ED 实现中
使用相同的占据截断进行对比。测试还验证了将 `Delta` 变为 `-Delta` 的显式组分相位旋转，
以及从化学势和逆温度的数值导数中得到的粒子数和巨正则能量。

高斯检查在 `beta=1`、`mu=-5`、无相互作用的周期 3×3 三角晶格上
运行六次独立计算，每次包含四个测量 bin。它们覆盖 `Delta=0,+0.2,-0.2`，
包括两个求解器之间的零配对对比。
每个 bin 的密度、物理能量和两个 K 点结构因子都与解析期望值进行比较；
配对求解器还测量了反常振幅。
这些测试同时验证了密度结构因子中的玻色恒等项和反常缩并项。
对比中不引入任何 benchmark 测量值或 ED 占据截断。

固定场检查开启两个相互作用通道（`U1=0.1`、`U2=-0.7`），
并在六个虚时片上读入显式非均匀场。
独立的稠密传播给出行列式权重和 Wick 缩并后的可观测量。
零提议位移固定辅助场不变，让可执行程序沿两个时间方向遍历，
取 `Nwrap=1,2,6`。程序还会检查最终场与给定场配置一致。
三次独立计算耗时约 2 秒。

改变 `Nwrap` 可检查二次型和给定有相互作用场算例中稳定化区间的一致性。
这些检查针对固定场传播和权重；
接受率决策以及严重病态低温乘积需要在相应参数下进行专项测试。

六次高斯计算在参考工作站上使用已编译的 Intel MPI/MKL 可执行文件
耗时约 11 秒。若只运行此检查（无需 QuSpin）：

```bash
BAFQMC_RUN_MPI_TESTS=1 python3 -m unittest discover -s tests/physics -p test_gaussian_solvers.py -v
```

不设置环境变量时，直接命令会跳过实机计算；`make physics` 会开启这些检查。
报告时请读取实际的测试摘要，确认哪些检查真正运行了。

## 采样和算法修改

高斯和有限 Fock 对比是确定性的，其容差反映浮点运算和有限差分精度。
有相互作用的蒙特卡洛计算同时具有采样误差以及 Trotter 误差和参考截断效应，
需要在所选物理参数下单独评估。

对于 HS 场、局域更新或有相互作用传播的修改，需要用相关可观测量和独立参考值
运行一次小型有相互作用计算。两个求解器都提供了实机回归计算任务：

```bash
make -C src/number_conserving benchmark-dqmc
make -C src/pairing benchmark-dqmc
```

粒子数守恒计算任务包含四个 100000-bin 的有相互作用算例，
在参考工作站上耗时约 10–15 分钟。
配对计算任务检验有限配对和零配对极限。
若修改了相互作用、几何结构或低温参数区间，需要补充相应的测试算例。

报告时需同时给出均值、分块 SEM、参考值和残差，以及采样和截断参数。
论文复现保留所有有效数据及其不确定度。
对于精确的可观测量约束使用确定性恒等式，对于采样对比使用实际不确定度进行解读。

## 工作流修改和新测试

对生产执行或续跑逻辑的修改还需要运行同时包含 BAFQMC 和 ED 的小型实机流水线测试：

```bash
BAFQMC_RUN_MPI_TESTS=1 python3 -m unittest discover -s tests/reproduction -p test_production.py -v
```

若 ED 使用独立环境，请设置 `BAFQMC_PYTHON_ED=/path/to/quspin/python`。
新的物理测试应标明算符或极限情形，计算独立期望值，并测验可能出错的数值输出。
显式保留种子、采样量、归一化以及任何截断参数。
生成的测量数据放入临时或被忽略的输出目录；
紧凑型测试输入和参考期望值可以提交到仓库中。

完整的论文复现是独立的生产计算，通过 `python3 reproduce.py` 调用。
其[资源指南](./benchmarks-paper-resources.md)给出了完整的计算预算。
