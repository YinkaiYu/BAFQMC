# 版本发布

每个版本对应一个计算源代码快照、其文档化的 benchmark 输入以及验证结果。
使用带版本号的 Git 标签和 GitHub Release，以便研究人员日后回溯到相同的代码状态。

## 准备版本

保持安装指南、Agent 指令和数值资源估计与代码的一致性。
检查源代码包是否包含 MIT 协议、求解器源码、ED 代码、精选输入、
处理后数据以及复现脚本。原始链数据和生成的编译文件不纳入 Release。

在已记录的 Linux 环境中运行以下检查：

```bash
python3 scripts/doctor.py
make check
make physics
make smoke
python3 reproduce.py --mode check
git diff --check
```

GNU Fortran 的 GitHub Actions 检查也应在 Release 提交上通过。
若生产编排逻辑有变更，还需运行真实的续跑测试：

```bash
BAFQMC_RUN_MPI_TESTS=1 python3 -m unittest discover -s tests/reproduction -v
```

对于物理修改，将哈密顿量、估计量归一化、Trotter 步长、采样统计量和
ED 占据截断与上一版本进行比较。运行[开发指南](./development.md)中描述的
相应小体系恒等式和实机回归计算。
除非对模型、输入设置或参考计算进行了明确的文档化更新，
否则应保持已发布的 benchmark 输入和处理后数据完整不变。

默认的 `python3 reproduce.py` 命令计算全部 22 个生产参数点，
预算约 **12–24 小时、16 GiB 内存和 8 GiB 可用磁盘**（参考桌面）；
详见[实测资源指南](./benchmarks-paper-resources.md)。
若版本修改了生产物理内容或 benchmark 输入，需重新运行受影响的生产扫描，
并记录其参数和结果。在 Release 说明中区分已完成的检查和完整计算。

## 发布源代码快照

提交所有变更后检查 `git status --short`。使用新版本号标记经过验证的确切提交，例如：

```bash
release_tag=v0.1.3
git tag -a "$release_tag" -m "BAFQMC ${release_tag}"
git push origin main
git push origin "$release_tag"
```

Source Release 工作流会对该标签重新运行计算 CI，
然后以 GitHub Release 的形式发布源代码存档。
在 Release 说明中包含 changelog 条目、相关验证结果，
以及对输入/输出格式的任何变更说明。标签一经发布应保持不变；
如需修正，发布新版本。

GitHub 会为带标签的快照提供源代码存档。
在将 Release 视为完成之前，请将存档下载到新目录，
并检查 `python3 reproduce.py --plan` 和
`python3 reproduce.py --mode check` 的运行情况。
提取出的包应在不依赖手稿检出或同级仓库的情况下正常工作。
可选的开发容器提供相同的 GNU/MPI 环境；
打开容器后会自动运行环境诊断，生产计算则由研究人员自行控制。
