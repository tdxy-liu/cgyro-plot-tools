# 压缩 FULL_T_ASYM 的跨机器读取（v0.2.20）

## 更新已经能启动的 Comparison 工具

在工具仓库目录中，用原来启动 GUI 的同一个 Python 环境运行：

```sh
git pull --ff-only
python -m pip install -r requirements-compression.txt
python cgyro_runtime_check.py
python cgyro_comparison.py
```

Windows 与 Linux 使用同样命令。若 Git 提示本机有未提交改动，请先检查、保留自己的改动，不要使用强制覆盖。
也可使用工具内的更新入口；更新源码后仍需在该 Python 环境安装可选的 zstandard 依赖并重启。

工具现在随仓库提供 `cgyro_fullt_reader.py`：
- FULLT 原格式、FTZ v1 压缩图、时间窗口、跨 Kx 追踪和文本导出统一使用它。
- 不再要求外部 pygacode 提供定制的 `fullt_compressed.py`，也不使用某台机器的绝对路径。
- 普通场、通量、其他诊断仍使用 pygacode；本次不替换其安装，不修改 CGYRO 求解器。
- 压缩解码需要 zstandard；只读取原始 FULLT 的路径不会导入它。
- 不提供假数据降级：模块缺失、数据损坏、格式不符均明确报错。

### 新机器尚未配置 pygacode

先安装原有 GUI 依赖：Python 3.9+、numpy、scipy、matplotlib、tkinter，以及 GACODE 的 pygacode。
本仓库不包含整个 GACODE。通常原来的绘图环境已经具备这些依赖。

可用 `GACODE_ROOT` 指向该机器的 GACODE 根目录（内含 `f2py/pygacode/cgyro/data.py`）。
也可在本工具目录创建 `cgyro_runtime.local.json`：

```json
{"gacode_root": "../gacode"}
```

相对路径从工具目录解析，而不是终端当前目录。
显式 `GACODE_ROOT` 优先，其次本机 JSON，最后是原有相邻目录/已安装的 pygacode。
配置错误时不自动换用另一份；更换后重启 GUI。
不需要把所有机器都指向“带压缩模块的新版 GACODE”；压缩模块已经随本工具分发。

## 离线安装

zstandard 是随 Python/操作系统编译的扩展，不能复制 Windows 的 .pyd 到 Linux，也不能混用 Python 3.9 和 3.12 的包。
`cgyro_runtime.local.json` 和 `.runtime/` 不上传 GitHub。

先在目标机器确认环境：

```sh
python -c "import sys,platform; print(sys.version); print(platform.machine()); print(platform.libc_ver())"
python -m pip --version
```

若目标是 **CPython 3.9 / Linux x86_64 / glibc 2.17**，在可联网机器下载目标平台的 wheel：

```sh
python -m pip download --only-binary=:all: --no-deps --implementation cp --python-version 39 --abi cp39 --platform manylinux2014_x86_64 --dest wheelhouse "zstandard==0.25.0"
```

把 wheelhouse 传到目标机器，然后在目标 CGYRO Python 环境内运行：

```sh
python -m pip install --no-index --find-links ./wheelhouse -r requirements-compression.txt
python cgyro_runtime_check.py
```

其他 Python/平台需替换下载参数。若软件源提示 versions: none，不能据此认定 CPU 或 Python 不兼容：也应检查服务器的软件源/联网限制。
没有合适 wheel 时需在目标平台构建，不要强装错误架构包。
高级用法仍支持 `.runtime/python/<cache_tag>-<platform>/` 私有依赖目录，只加载与当前解释器匹配的目录；常规安装不需要此目录。

## 算例文件与选择规则

压缩算例保留：
- `bin.cgyro.fullt_asym.ftz` 与同一文件身份的 `.fti`；
- `input.cgyro` / `input.cgyro.gen`、网格、时间标签及普通 GUI 所需输出。

`FULL_T_ASYM_COMPRESSION=1` 选择 FTZ，0 选择原始文件。
两者并存时不会混读；转换后遗留的原始副本不会随压缩续算更新。
数据位于算例根目录或 bin 子目录均支持，但同格式两处同时存在会拒绝歧义选择。
选择与实际文件冲突、数据/索引身份不符、布局或时间标签不符时明确报错。

原始实数 FULLT 的 MPI 布局按运行时 `TOROIDALS_PER_PROC` 还原，不按文件大小猜 reshape。
需要运行对应的原输入、网格和完整时间标签；旧复数容器不能静默按新实数布局读入。
当 N_RADIAL 不能整除 BOX_SIZE 时，按原生整数除法读取 thetab 长度并保留分辨率警告，不改变网格。

## 按需读取、缓存与导出

单 K 图仅解压所选被诊断 Ky 和时间的块；每块包含全部已保存 Kx，因此选择单 Kx 也需在内存解码对应块，但不会写出完整未压缩文件。
追踪只保存请求时间窗的跨 Kx 数据。大时间窗仍可能占用较多内存。
FTZ 仅读取完整提交的索引记录；追加、回退或文件替换会改变缓存签名。已提交损坏不作为“未完成尾部”忽略。
实时输出时索引与时间文本数量可短暂不同，但重叠时间前缀必须一致。原格式没有提交索引，文件尺寸不完整会报错。

导出列表把 .ftz 识别成 FULL_T_ASYM，忽略 .fti 和 partial 文件。
文本按块导出到用户指定的新文件，不整文件落盘解压；全文本可能远大于二进制，不建议磁盘紧张时全量导出。
已有的文本路径不会由该导出器自动覆盖。

```sh
# 检查环境、索引和一个已提交数据块（不是全量数据验收）
python cgyro_runtime_check.py --case /path/to/case

# 显式解码校验全部块（耗时取决于数据量；原始 raw 没有内置校验和）
python cgyro_runtime_check.py --case /path/to/case --verify-all
```

## 科学含义未改变

此次更新只涉及读取/分发，不改符号、精度、归一化、阈值、色标或诊断频率。
GUI 中 source 是历史接口命名；固定 K 是“被诊断模”，不是唯一源模。
无损解码不强制配对对称化，不删除边界，不补造数据。
现有逐图色标和显示阈值仍需在物理比较中独立审查，不能用颜色直接比较绝对强度。
原诊断的 RHS 采样时刻和径向 Nyquist 别名等约定不会因读取更新而改变。

## 测试与范围

```sh
python -m unittest discover -s tests -v
```

普通测试不依赖外部 pygacode 的压缩模块；3 项定制 GACODE 集成测试默认跳过。
拥有对应的更新版 GACODE 时，可设置 `CGYRO_TEST_GACODE_INTEGRATION=1` 额外测试 pygacode 全量提取和两个 GACODE 验证器。
GUI 的正常读取不通过这些外部压缩接口。

逐位测试覆盖 float32/64、单/双通道、Kx0/全 Kx、正负零、极小数、NaN/Inf、非对称配对、原样块、损坏提交、未提交尾部及回退。
原格式切片/追踪/导出覆盖 TOROIDALS_PER_PROC=1/2/4。
这些是后处理测试，不等同于再次通过 Linux CGYRO 求解器/MPI/自适应时间步的运行验收。
读取器来源和许可证见 `third_party/gacode/`。
