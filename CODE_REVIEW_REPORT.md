## Code Review Report

审查对象：CGYRO Comparison Tool v0.2.22 源码，基线提交为 fa14080；本轮修复发布版本为 v0.2.23。

审查日期：2026-10-29。

范围包括 GUI、case 与 workspace 管理、通量与涨落绘图、FFT、能量与 FULLT 读取、数据导出、更新与代理、远程提交脚本以及测试。下列 Critical 和 Warning 项已修复，并纳入 v0.2.23 提交。用户原有的 USER_GUIDE.md 和 USER_GUIDE_legacy.md 改动保留在本地。

### Critical

#### C1 运行记录可执行任意 Shell 命令

- 位置：remote_cgyro 的 load_run_info、SSH 远端命令和 rsync 调用。
- 问题与原因：原先使用 `. "$RUN_INFO_FILE"` 读取运行记录；修改或分享的 `.env` 文件可以在拉取结果前执行本地命令。部分远端命令还把目录直接拼进命令字符串，空格、引号、分号会改变命令含义。
- 修复：只接受运行记录的固定字段，按 printf %q 的反斜杠格式解析值；不再 source/eval 记录。检查记录地域，使用 printf %q 引用 SSH 路径，并让 rsync 通过 --protect-args 传递路径。
- 兼容范围：保留脚本生成的普通路径、空格和转义引号；不接受控制字符路径或手写的 Shell 表达式。
- 修复示例，字段白名单验证之后执行：

```bash
decoded=$(decode_run_info_value "$encoded") || return 1
printf -v "$key" '%s' "$decoded"
```

- 验证：隔离测试确认普通命令与命令替换不会执行，旧格式路径可以恢复，地域不匹配会拒绝；测试未连接 SSH、rsync 或 Slurm。

#### C2 workspace 保存失败会破坏原文件

- 位置：cgyro_comparison_ui.py 的 _write_workspace_file；cgyro_update.py 的 save_update_proxy_config。
- 问题与原因：workspace 以 w 模式直接打开，序列化失败或写入中断时旧文件已被截断。代理设置使用固定名称的 .tmp 文件，也存在并发写入相互覆盖的风险。
- 修复：新增 cgyro_io.py，先完成 JSON 序列化，再写同目录的唯一临时文件，flush/fsync 成功后用 os.replace 替换。失败时清理临时文件并保留旧文件。当前图文本导出也复用此写入方式。
- 修复示例：

```python
content = json.dumps(payload, allow_nan=False)
with atomic_text_file(path) as handle:
    handle.write(content)
```

- 验证：不可序列化值和模拟替换失败均不会修改旧 workspace，不遗留临时文件。

#### C3 加载失败的 workspace 会清空当前工作状态

- 位置：cgyro_comparison_ui.py 的 _load_workspace_file、_restore_auto_workspace 和更新重启失败处理。
- 问题与原因：旧流程先清空 case 列表，再尝试加载新数据。共享盘断开、路径不可用或 selected_cases 格式错误时，当前 case、选择和 Advanced 状态可能先被丢弃。保存的 case 别名也没有被恢复。
- 修复：先验证结构、解析相对路径并暂存新 case，再替换当前 workspace；全部读取失败时保留当前状态。恢复 case 别名及选择，部分失败时报告跳过项；自动恢复或重启失败时保留恢复文件。
- 修复示例：

```python
if cases and not staged:
    return None
self.cases = staged
```

- 验证：全部失败、部分失败、无效选择结构、同名目录的不同别名、相对路径、自动恢复失败均有测试。

#### C4 二进制可读导出可能覆盖已有结果

- 位置：cgyro_data_export.py 的 transfer_bin_to_readable 和 _new_readable_export_dir。
- 问题与原因：重复使用 case_readable 目录，并以 w 模式写文件。重复导出会覆盖上次结果，不同名称清洗后相同也会共用目录。
- 修复：以原子 mkdir 创建新目录，已存在时依次分配 _2、_3 等后缀；单个 case 失败时继续处理其他 case 并报告错误。
- 修复示例：

```python
try:
    os.mkdir(path)
    return path
except FileExistsError:
    index += 1
```

- 验证：重复导出和名称清洗冲突分配不同目录，旧结果保持不变。

### Warning

#### W1 缺失精度属性时误用 float64

- 位置：cgyro_data_export.py 的 _real_dtype_for_data。
- 原因：np.dtype(None) 会返回 float64，因此缺失 BYTE 属性反而提前选中了双精度，未进入所注释的 float32 后备路径。直接读取单精度文件时，值和元素数都会错误。
- 修复：显式跳过 None，依次使用当前 case 的精度标记、实例或模块 BYTE，最后才使用 float32；同时支持多列 equilibrium 文件。
- 示例：`if dtype is None: continue`。
- 验证：缺失属性、HIPREC、当前 case 标记优先级及实际 float32 文件的数值均有测试。

#### W2 掩码丢失和图像上下方向错误会改变导出数据

- 位置：cgyro_data_export.py 的 _record_current_plot_xyz_dataset、_flatten_xyz_dataset 和 _collect_current_plot_xyz_datasets。
- 原因：先 np.asarray 再 filled 会丢掉 MaskedArray 的掩码；图像回退导出也没有区分 origin='upper'，可能把数值配给错误的 y 坐标。
- 修复：用 np.ma.asarray 保留掩码后填充 NaN；按图像 origin 排列 y 坐标。
- 示例：`z = np.ma.asarray(z, dtype=float).filled(np.nan)`。
- 验证：掩码二维图、掩码曲线、upper 图像的坐标和值对应关系。

#### W3 当前图导出包含重复或不可见数据

- 位置：cgyro_data_export.py 的 _collect_current_plot_xyz_datasets；FFT 与 Flux 的显式缓存调用。
- 原因：缓存过的 imshow 数据还会被艺术对象回退逻辑再次收集；隐藏曲线和 colorbar 也未被过滤。
- 修复：缓存记录所属 axes，缓存图像不再二次导出，跳过隐藏对象及 colorbar；多面板 FFT/Flux 明确传入所属 axes。
- 示例：`if ax in cached_axes or not img.get_visible(): continue`。
- 验证：缓存图像只导出一次，隐藏曲线和 colorbar 不进入数据集。

#### W4 导出元数据可能与当前图不一致

- 位置：cgyro_data_export.py 的 _capture_current_plot_context、_current_plot_export_metadata；Plotting.plot_comparison。
- 原因：绘图后修改时间窗口或 case 选择，再直接导出，旧实现会用新控件状态解释旧图。
- 修复：绘图时深拷贝模式、case 和 workspace 状态；导出沿用绘图快照。坐标轴标题仍取当前图本身。
- 示例：`metadata.update(copy.deepcopy(self._current_plot_context))`。
- 验证：绘图后把 Time Start 从 500 改为 900，导出图的元数据仍记录 500。

#### W5 旧 Release 会挡住主分支的新版本

- 位置：cgyro_update.py 的 check_for_updates。
- 原因：只要 Release API 返回合法版本就立即返回；例如当前为 0.2.22、Release 为 0.2.4、main 为 0.2.23 时，会误报没有更新。
- 修复：比较 Release 与主分支 VERSION，选择较新版本；单个服务失败时使用另一个有效结果。
- 示例：`max(candidates, key=lambda info: _version_key(info.latest_version))`。
- 验证：旧 Release、新 Release、单服务离线、两个服务均失败等路径。

#### W6 Git 更新阻塞 GUI 且后台线程访问 Tk

- 位置：cgyro_comparison_ui.py 的 _pull_update_from_git、_post_update_check_result 和新增的 worker/poll 方法。
- 原因：Git 下载在 Tk 主线程执行，可能冻结窗口至超时；检查线程直接 root.after 也会在窗口退出或部分 Tk 环境中触发跨线程错误。
- 修复：网络与 Git 操作运行于后台线程，只向 queue.Queue 写结果；Tk 主线程轮询队列，更新控件或重启。失败统一释放忙状态。
- 示例：`self._update_events.put((callback, args))`，由 _poll_update_events 调用 callback。
- 验证：后台检查不调用任何 Tk API；阻塞模拟下载期间启动方法能返回；意外异常也能释放更新状态。

#### W7 不规则采样被当作等间距 FFT

- 位置：cgyro_comparison_plotting_fft.py 的 _prepare_fft_samples 和 _plot_fluctuation_fft。
- 原因：旧实现使用时间差绝对值的平均值，倒序、重复或不规则采样仍可能得到看似正常的频率轴；数据输出较短时也可能按时间索引越界。
- 修复：先对齐数据与时间长度，只使用有效索引，要求有限、严格递增且近似等间距的采样；不满足时明确拒绝 FFT，提示先重采样。
- 示例：`field, times, dt = self._prepare_fft_samples(field, data.t, t_indices)`。
- 验证：不同长度、非均匀、重复、倒序和 NaN 时间轴。

#### W8 固定阈值使小振幅数据无法归一化

- 位置：cgyro_comparison_plotting_fluctuation.py 的 _normalize_fluctuation_profile 和自定义分母；UI 自定义归一化输入。
- 原因：把绝对值小于 1e-15 的非零分母当成零，1e-20 量级的线性涨落不能按 max/min 归一化，且没有对应归一化标签。
- 修复：检查分母是否有限和是否真正为零，不使用带绝对量纲的固定阈值。
- 示例：`if not np.isfinite(denominator) or denominator == 0.0: ...`。
- 验证：1e-20 的 max/min/min-max，1e-30 的自定义分母，以及全零数据。

#### W9 重复把整个复数场重新转换

- 位置：cgyro_comparison_plotting.py 的 _coerce_kxky_complex 和 _load_kxky_complex。
- 原因：每次绘图都执行 real + 1j*imag，重复分配完整的四维或五维场。若直接替换 pygacode 公共数组，又会破坏 kxky_select 对 packed 布局的约定。
- 修复：标准 Fortran 数据采用共享内存的复数视图；其他布局缓存转换结果。以原数组身份失效缓存，保留 pygacode 公共数组。
- 示例：`np.ndarray(arr.shape[1:], dtype=dtype, buffer=arr, strides=arr.strides[1:])`。
- 验证：重复加载返回同一转换结果，Fortran 视图共享内存，原数组不变，更换源数组后重新转换。

#### W10 径向相关图不必要地读取全部大场

- 位置：cgyro_comparison_plotting_other.py 的 _plot_other_rcorr_phi。
- 原因：仅需一个 Phi/Apar/Bpar 切片，却可能调用 getbigfield 加载所有矩；先画 Phi 再画 Apar/Bpar 时，单纯检查 kxky_phi 是否存在也不能保证目标场已经加载。
- 修复：经共享场读取器只加载选中的场，再在内存中选 theta，不再依赖 getbigfield/kxky_select 的全局状态。
- 示例：`field = self._load_named_kxky_complex(data, label, field_name)`。
- 验证：选择 Apar 时只读取 Apar，两个全量接口均不调用。

#### W11 无变化的选项刷新重复绘图并读取 case 元数据

- 位置：cgyro_comparison_ui.py 的 _draw_formula_panel；cgyro_comparison_plotting.py 的 _cached_input_scalars 和 _read_input_cgyro_scalars。
- 原因：选项事件会反复清空并绘制相同公式；扫描选项还会重新读取 input.cgyro，共享文件系统上会放大延迟。
- 修复：按公式内容、DPI 和显示参数复用画布；标量读取使用有界缓存，文件时间与大小变化即失效，对调用者返回字典副本。
- 示例：`if getattr(canvas, '_formula_render_key', None) == render_key: return`。
- 验证：相同内容连续刷新三次只渲染一次，内容变化再刷新；重复读文件只打开一次，文件改动会重新读取。

#### W12 更新连接参数可能被当作命令选项

- 位置：cgyro_update.py 的代理解析、_validate_update_proxy_config、SshSocksTunnel.start 和 update_from_git。
- 原因：argv 列表本身不能阻止 SSH 把以 '-' 开头的主机解释为选项；SOCKS URL 也接受端口 0、路径和部分不适合进入代理子命令的字符。
- 修复：验证主机、用户名、远端和分支参数；SSH 目标前加入 '--'；拒绝非法 SOCKS URL，规范 IPv6 URL。
- 示例：`if remote.startswith('-'): raise GitUpdateError(...)`。
- 验证：SSH/Git 选项输入在进程启动前拒绝，非法 URL、Shell 元字符和 IPv6 往返均有测试；脏目录、错误分支和非 fast-forward 仍会拒绝更新。

#### W13 命令行更新依赖 GUI 和科学计算环境

- 位置：cgyro_comparison.py 入口；FULLT bridge、energy 和 export 的导入。
- 原因：原入口先导入 Tk、NumPy 和 pygacode，再解析 --update；环境缺失时，连更新入口也不能启动。部分 FULLT 路径使用绝对导入，破坏已支持的 package 导入方式。
- 修复：先分流命令行更新参数，再导入 GUI；FULLT 导入统一支持脚本和包模式。
- 示例：在 GUI imports 之前处理 `if __name__ == '__main__' and len(sys.argv) > 1`。
- 验证：python -S cgyro_comparison.py --help 不加载第三方依赖即可运行；从父目录导入 CGYRO.cgyro_comparison 成功。

#### W14 CSV 标签可能被表格软件当作公式

- 位置：cgyro_data_export.py 的批量 CSV 导出和 _spreadsheet_safe_label。
- 原因：case/曲线标签是文本输入，开头的 =、+、-、@ 在表格软件中可能被当作公式；CSV 引号本身不能把它变为普通文本。
- 修复：对这类 Dataset 标签加文本前缀；数值列和 JSON 元数据中的原始标签保持原义。
- 示例：`"'" + text if text.lstrip().startswith(('=', '+', '-', '@')) else text`。
- 验证：公式样式标签转成文本，普通 case 标签不变。

### Info

#### I1 通用职责已抽取，能量模块仍适合分阶段重构

- 位置：cgyro_io.py、GUI 更新 worker、FFT 前处理，以及 cgyro_comparison_plotting_energy.py。
- 本轮改进：统一原子文件写入；分开网络操作和 UI 事件处理；提取可独立测试的 FFT 采样检查及缓存入口，减少隐式依赖。
- 后续可选：基线约有 19,715 行 Python 源码，38 个方法超过 100 行；FULLT 主加载方法约 691 行，绘图方法也较长。建议在真实算例的数值与导出回归基准下，继续按读取、坐标变换、归一化和绘制拆分。这些长函数没有仅为行数目标而被机械改写。
- 后续接口示例：`diagnostic = load_diagnostic(case, selection)`、`profile = reduce_diagnostic(diagnostic, window)`、`render_profile(axis, profile)`。

#### I2 GUI 检查原先在常规测试中长期跳过

- 位置：.github/workflows/fullt-portability.yml、tests/test_review_regressions.py 和 tests/test_remote_launcher.py。
- 改进：新增离线、无数据依赖的故障回归与隔离 Bash 测试；现有 CI 增加 Linux Xvfb 和 Windows 的 Tk 布局步骤。
- 本机已运行隐藏 Tk 窗口测试；新增的远端 CI 步骤需在提交后由 GitHub Actions 执行。

### Summary

- Critical：4 项，已修复。
- Warning：14 项，已修复。
- Info：2 项；通用职责和 CI 已改善，较长能量/FULLT 方法保留为后续可选重构。
- 总体评价：主要风险集中在状态保存、命令边界、二进制精度和导出正确性。已有的模块划分与 FULLT 完整性测试有价值；本轮增加了错误路径、边界输入和 GUI 线程行为的覆盖。

#### 六个审查维度

| 维度 | 主要结论 |
| --- | --- |
| 可读性与风格 | 新增 helper 采用 snake_case，解释缓存失效、原子替换与 Fortran 视图的原因；仍有较长能量方法 |
| 架构原则 | 共用文件写入，分离 worker/UI，保留 pygacode 原始布局；GUI 和数值后端仍通过较多 Mixin 属性协作 |
| 测试覆盖 | 新增 46 项测试，覆盖保存失败、加载失败、代理与参数验证、导出、FFT、归一化、缓存和 Shell 记录解析；未测得行/分支覆盖率百分比 |
| 安全性 | 修复可执行运行记录、远端路径引用、命令选项混淆和 CSV 公式问题；未发现本项目硬编码密码或 Token |
| 性能 | Git 不阻塞 Tk；减少公式重绘、input 元数据重读和整场复数转换；径向相关只读取目标场 |
| 可维护性 | 失败保留 workspace/恢复文件，导出避免覆盖，导出元数据绑定绘图状态；CI 增加布局检查 |

#### 验证方法

本机使用 Windows Python 3.12，运行完整 unittest 和隐藏 Tk 布局测试：

```powershell
$env:CGYRO_TEST_TK_LAYOUT = '1'
python -m unittest discover -s tests
```

最终测试：96 项，93 项通过，3 项按配置跳过。跳过项需要更新版外部 GACODE，未冒充本机已经验证。另检查了 Python 源码语法、Bash 语法和本轮代码的 diff 空白。

现有测试包含合成数值、native raw/FTZ 读写、CRC/损坏数据、时间追加和缓存失效。未使用用户的真实大体积算例做全量数值重算，也未执行真实集群提交或 GitHub 更新下载。远端 rsync 需支持 --protect-args。

依照 [OSV 批量查询接口](https://google.github.io/osv.dev/api/)，本机查询了 numpy 2.1.1、matplotlib 3.9.2、scipy 1.14.1、zstandard 0.25.0 四个精确版本，响应均未返回已知漏洞记录。该结果限于这四个 Python 包，不代表自定义 pygacode、系统原生库和所有传递依赖已完成漏洞审计。

本轮文档沿用项目内 Markdown 格式，按 write-page 的写作和核验规范整理。修复随 v0.2.23 源码提交；GitHub Release 单独管理。
