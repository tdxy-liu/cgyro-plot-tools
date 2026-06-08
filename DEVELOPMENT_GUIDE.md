# CGYRO Comparison Tool — 超详细开发说明书

## 1. 项目总览

**CGYRO Comparison Tool** 是一个基于 `Tkinter` + `Matplotlib` 的多算例交互式对比可视化工具。核心任务是：加载一个或多个 CGYRO 输出目录，解析 `pygacode` 数据结构，并在同一画布上对比绘制频谱、通量、涨落、FFT、带状流（Zonal ExB）、能量平衡、径向相关、POD 奇偶分解等诊断量。

**设计原则**：
- **Mixin 组合模式**：UI、数据导出、绘图基类、各专科绘图全部解耦到独立文件。
- **防御式编程**：`pygacode` 导入失败时降级为 Mock；所有数组操作先做 `size/shape/ndim/finite` 校验。
- **时间窗统一抽象**：所有需要时间平均的绘图共享 `_get_time_indices(t)`。
- **全局视觉风格统一**：颜色板、字体、线宽全部集中在 `cgyro_comparison_plotting.py` 顶部常量。

---

## 2. 完整文件结构与模块地图

| # | 文件 | 行数 | 核心职责 | 导出符号 |
|---|------|------|----------|----------|
| 1 | `cgyro_comparison.py` | 29 | **程序入口**。继承链终点 `CGYRO_Comparison`；创建 `tk.Tk()` 并进入主循环。 | `CGYRO_Comparison` |
| 2 | `cgyro_comparison_bootstrap.py` | ~80 | **引导层**。`pygacode` 导入适配；失败时提供 `Mock` 降级；全局常量与默认目录 helper。 | `cgyrodata`, `cgyrodata_plot`, `default_share_dir`, 常量 |
| 3 | `cgyro_comparison_ui.py` | ~2600 | **UI 与状态机**。窗口布局、菜单、集中式下拉选项常量、所有控件（含 40+ 个动态选项控件）、用例增删改查与拖拽排序、时间窗控制、公式面板渲染系统、物种解析系统。 | `CgyroUiMixin` |
| 4 | `cgyro_comparison_plotting.py` | ~2518 | **绘图调度中枢**。全局样式、时间窗解析、1D 绘图包装、输入文件 diff、用例信息滚动浏览、POD/ZF/连通kx 数学工具、 `_plot_single_case` 路由表。 | `Plotting` |
| 5 | `cgyro_comparison_plotting_frequency.py` | ~78 | **频率 / 增长率**。`freq` 数组时均与 `/ky` 归一化。 | `FrequencyPlotting` |
| 6 | `cgyro_comparison_plotting_flux.py` | ~801 | **通量诊断**。`ky_flux` 的 vs ky / vs Time / vs ky_time / vs kx(估计) / vs 2D 扫描、场分解、实离子 GyroBohm 归一化。 | `FluxPlotting` |
| 7 | `cgyro_comparison_plotting_fluctuation.py` | ~559 | **涨落诊断**。1D（vs ky/kx/Time/fft）、2D 实空间（vs xy 动画 / vs xt 等高线）、`vs kxky` 频谱平面图、电子尺度空间归一化。 | `FluctuationPlotting` |
| 8 | `cgyro_comparison_plotting_fft.py` | ~372 | **FFT 谱**。沿时间轴 FFT，支持 Linear coherent / Nonlinear incoherent、Amplitude / Power、线性频率叠加。 | `FftPlotting` |
| 9 | `cgyro_comparison_plotting_zf.py` | ~410 | **带状流**。vs Time / vs kx / phi vs kx(theta=0) / vs gamma_lin（含比值计算）。 | `ZfPlotting` |
| 10 | `cgyro_comparison_plotting_energy.py` | ~1679 | **能量平衡**。基于 `triad_v2` 的熵平衡、ZF 能量平衡、有效增长率 gamma_eff、单图、FULLT 传递图、2D 扫描。 | `EnergyPlotting` |
| 11 | `cgyro_comparison_plotting_other.py` | ~523 | **其他诊断**。积分误差、rcorr_phi（含指数包络拟合）、POD 奇偶分解（SVD + 镜面对称分析）。 | `OtherPlotting` |
| 12 | `cgyro_data_export.py` | ~332 | **数据导出**。`bin.cgyro.*` 批量转制表符文本；`.cgyro.triad` 按 8 通道拆 Sheet。 | `CgyroDataExportMixin` |

---

## 3. 类继承与 Mixin 体系（精确签名）

```python
class CGYRO_Comparison(CgyroUiMixin, CgyroDataExportMixin, Plotting):
    pass
```

### 3.1 CgyroUiMixin（状态机 + 视图层）

**状态属性（在 `__init__` 中初始化）**：

| 属性 | 类型 | 说明 |
|------|------|------|
| `self.cases` | `dict[str, cgyrodata]` | 已加载算例的中央仓库，key 为 Listbox 显示名 |
| `self.ani` | `FuncAnimation \| None` | 当前 Matplotlib 动画对象 |
| `self.is_paused` | `bool` | 动画暂停标志 |
| `self.current_frame` / `self.total_frames` | `int` | 动画/手动翻页状态 |
| `self._drag_start_index` | `int \| None` | 算例列表拖拽排序的起始索引 |
| `self._manual_pager_active` | `bool` | 是否在手动翻页模式 |
| `self._energy_entropy_axes_active` | `bool` | 能量平衡 entropy vs ky 是否启用了双轴模式 |
| `self._energy_entropy_axes` | `tuple[Axes, Axes] \| None` | 双轴缓存（上=ion，下=electron） |
| `self.fluc_average_var` | `StringVar` | 平均模式：`"Root Mean Square"` / `"Mean Absolute"` |
| `self.time_mode_var` | `StringVar` | 时间菜单模式：`"Manual Start/End"` / `"Last percent"` / `"Last duration"` / `"Full range"` |
| `self.time_percent_var` | `StringVar` | `"50"`，Last percent 模式使用 |
| `self.time_duration_var` | `StringVar` | `""`，Last duration 模式使用 |
| `self.plot_type_var` | `StringVar` | 当前 Plot Type 组合框值 |
| `self.t_start_var` / `self.t_end_var` | `StringVar` | 左侧时间输入框 |
| `self.log_x_var` / `self.log_y_var` | `BooleanVar` | 对数坐标开关 |
| `self.norm_ky_var` | `BooleanVar` | Frequency/Growth Rate 的 `/ky` 开关 |
| `self.flux_type_var` | `StringVar` | `"Energy"` / `"Particle"` |
| `self.flux_xaxis_var` | `StringVar` | `"v.s ky"` / `"v.s kx (estimated)"` / `"v.s Time"` / `"v.s ky_time"` / `"v.s 2D"` |
| `self.flux_decomp_var` | `BooleanVar` | 按场分解（Phi/Apar/Bpar） |
| `self.flux_norm_real_ion_var` | `BooleanVar` | 实离子 GyroBohm 归一化 |
| `self.fluc_field_var` | `StringVar` | `"Phi"` / `"Apar"` / `"Bpar"` |
| `self.fluc_xaxis_var` | `StringVar` | `"v.s ky"` / `"v.s kx"` / `"v.s Time"` / `"fft"` |
| `self.moment_var` | `StringVar` | Fluctuation 2D 的矩选择 |
| `self.fluc2d_view_var` | `StringVar` | `"vs xy"` / `"vs xt"` / `"vs kxky"` |
| `self.fluc2d_x_elec_var` | `BooleanVar` | 空间轴用电子尺度 `rho_e` |
| `self.zf_xaxis_var` | `StringVar` | `"vs Time"` / `"vs kx"` / `"phi vs kx(theta=0)"` / `"vs gamma_lin"` |
| `self.linear_gamma_file_var` | `StringVar` | 线性谱文件路径（ZF 和 FFT 共用） |
| `self.zf_gamma_lin_ky_var` | `StringVar` | vs gamma_lin 模式下用于插值比值的 ky |
| `self.energy_balance_mode_var` | `StringVar` | `"Entropy balance"` / `"ZF energy balance"` / `"Effective growth rate"` / `"Single plot"` / `"FULLT transfer map"` / `"2D scan"` |
| `self.energy_balance_n_var` | `StringVar` | n 索引（或 gamma_eff 模式下的 ky 值） |
| `self.energy_balance_spec_var` | `StringVar` | `"Total (-1)"` / `"Main ion (0)"` / `"Electron (1)"` |
| `self.energy_balance_single_quantity_var` | `StringVar` | `"T"` / `"N"` / `"T-N"` / `"entropy"` |
| `self.energy_balance_single_xaxis_var` | `StringVar` | `"vs Time"` / `"vs ky"` / `"vs kxky"` |
| `self.energy_balance_transfer_quantity_var` | `StringVar` | `"Re"` / `"Im"` / `"Abs"` |
| `self.energy_balance_transfer_ky_var` | `StringVar` | FULLT 模式固定 ky |
| `self.others_plot_var` | `StringVar` | `"Error"` / `"rcorr_phi"` / `"POD_parity"` |
| `self.others_rcorr_field_var` | `StringVar` | rcorr 的场选择 |
| `self.others_rcorr_theta_var` | `StringVar` | rcorr 的 theta 索引 |
| `self.others_pod_field_var` | `StringVar` | `"Apar"` / `"Phi"` |
| `self.others_pod_kx_var` / `self.others_pod_ky_var` | `StringVar` | POD 的 kx/ky 选择 |

**核心方法**：

- `_create_layout()` → 构建左右分栏主框架：左侧 `Canvas+Scrollbar` 内嵌控制面板；右侧 `FigureCanvasTkAgg` + `NavigationToolbar2Tk`。
- `_create_menu()` → 四级菜单：`File`（加载）、`Data`（导出）、`Time`（时间模式单选）、`Average`（平均模式单选）。
- `_init_options()` → **创建但不布局**所有动态选项控件（约 40+ 个 `ttk` 对象），并统一绑定 `<<ComboboxSelected>>` 到 `self.update_options`。这是关键：**控件只创建一次，通过 `grid_remove()` / `grid()` 切换可见性**，避免反复销毁带来的性能问题和事件泄漏。
- `update_options(event=None)` → 根据 `self.plot_type_var.get()` 决定显示哪些动态控件。使用 `row` 计数器从上往下逐行 `grid()`，隐藏其余。
- `_hide_dynamic_options()` → 对预定义的 `widgets` 列表（约 40 个控件引用）统一调用 `grid_remove()`。
- `_create_formula_panel(master, figsize, dpi)` → 在 `master` 内创建一个 **可滚动的 Matplotlib 小图**：`Canvas` → `create_window` 嵌入 `FigureCanvasTkAgg` 的 tk widget → 绑定 `Scrollbar`。返回 `(frame, fig, ax, canvas, viewport, vscroll, hscroll)`。
- `_draw_formula_panel(fig, ax, canvas, lines, widget=..., ...)` → **公式渲染引擎**。核心逻辑：
  1. 逐行分析文本，检测 `\frac`、 `\sqrt`、 `\sum`、 `\langle` 等 token，按行分配不同字号和行高（像素级）。
  2. 累加像素高度 → 计算所需 fig height → 调用 `fig.set_size_inches(fig_w, fig_h, forward=True)`。
  3. 用 `ax.text()` 逐行绘制在 `transAxes` 坐标系，`clip_on=False` 防止截断。
  4. 同步更新 Tk widget 的 `width/height` 和 Canvas 的 `scrollregion`，确保滚动条范围正确。
  5. 调用 `canvas.draw_idle()` 而非 `draw()` 减少闪烁。
- `_render_flux_kx_formula_math()` / `_render_fluctuation_1d_formula_math()` / `_render_fft_formula_math()` / `_render_zf_formula_math(mode_text)` / `_render_energy_balance_formula_math()` / `_render_pod_formula_math()` → 各模式调用 `_draw_formula_panel()` 的具体内容生成器，返回 LaTeX mathtext 字符串列表。
- `_safe_mathtext_line(line)` → 将 `\ge` 替换为 `\geq`，`\le` 替换为 `\leq`，兼容旧版 Matplotlib。
- `_formula_font_size(line, base_fs, ...)` → 根据视觉长度估算自动缩小字号，防止长公式在窄面板中被截断。

**算例管理方法**：

- `add_case_single()` → `filedialog.askdirectory` 单目录加载。
- `_select_multiple_case_dirs(initial_dir)` → 弹出独立 `Toplevel` 对话框，在父目录下列出子文件夹，支持 `Ctrl/Shift` 多选，返回选中目录列表。
- `add_case_multiple()` → 调用 `_select_multiple_case_dirs` 批量加载。
- `add_group()` → 选择父目录后自动加载所有 `_looks_like_case_dir(sub_path)` 的子目录。
- `_looks_like_case_dir(dir_path)` → 启发式判断：目录内是否存在 `input.cgyro`、`out.cgyro.freq`、`input.gacode` 任一文件。
- `_load_case_from_dir(dir_path, silent=False)` → 核心加载逻辑：
  1. 统一路径分隔符为 `/`。
  2. 生成唯一名称（重名时追加 `_1`, `_2`...）。
  3. 优先 `cgyrodata(dir_path)`，失败则降级 `cgyrodata_plot(dir_path)`。
  4. 存入 `self.cases` 并在 Listbox 追加条目。
- `remove_case()` / `remove_all_cases()` / `reload_cases()` → 删、清、重载。重载时保留原名称，从 `data.dir` 或 `data.path` 找回目录。
- `_on_drag_start(event)` / `_on_drag_motion(event)` → Listbox 内鼠标拖拽排序，直接操作 Listbox 的 `delete/insert`，并同步保持选中状态。

**物种解析系统**：

- `_get_case_species(data)` → 从 `data.z` 和 `data.mass` 提取 `(Z, M)` 列表。Mass 在 CGYRO 中已归一化到氘（~2.014 u）。
- `_get_species_name(z, m)` → 根据 Z 和 M*2（近似原子质量）识别化学元素名：Electron、Hydrogen、Deuterium、Tritium、Helium-3/4、Lithium...Tungsten，其余返回 `"Ion (Z=..., M=...)`。
- `_get_species_short_name(z, m)` → 返回紧凑符号：e、D、He3、He4、Li...W，未知则 `Z{int}`。
- `_get_case_species_densities(data, n_species)` → 从 `data.dens` 或 `input.cgyro` 的 `DENS_{i+1}` 获取归一化密度。
- `_get_main_ion_indices(data)` → 返回 `Z>0` 且 `DENS>0.3` 的物种索引。
- `_get_all_ion_indices(data)` → 返回所有 `Z>0` 的物种索引。
- `_update_species_list()` → 遍历所有已加载算例，取 **共同物种交集**，更新 `self.species_combo['values']`。优先添加 `"Main Ion (DENS>0.3)"` 和 `"All Ions"`。
- `_resolve_species_indices(data, n_species, case_label, species_override_index=None, fallback_first=True, main_ion_policy="all", single_species_only=False)` → **核心物种选择解析器**：
  - 如果 `species_override_index` 不为 None，直接返回该索引（用于 `"Plot All Species"` 循环）。
  - GUI 选择 `"Main Ion"` → 调用 `_get_main_ion_indices`；`main_ion_policy="first"` 时仅保留第一个主离子。
  - GUI 选择 `"All Ions"` → 返回所有离子索引。
  - GUI 选择具体物种（匹配 `"Z=..., M=..."` 文本）→ 精确匹配 Z/M。
  - 都失败且 `fallback_first=True` → 返回 `[0]`。
  - 最后受 `single_species_only` 约束：如果解析出多个索引，只保留第一个并打印提示。

**时间控制方法**：

- `_set_time_last_percent()` / `_set_time_last_duration()` → 弹出 `simpledialog.askfloat` 让用户输入百分比或物理时长，写入 `self.time_percent_var` / `self.time_duration_var`，并切换 `time_mode_var`。
- `clear_time_range()` → 清空所有时间相关变量，回到默认 `"Manual Start/End"`。
- `_clear_simple_time_entries()` → 仅清空 `t_start_var` / `t_end_var`，用于 `"Full range"` 模式。

**2D 扫描参数系统**：

- `_collect_flux_scan_case_scalars(selected_case_names)` → 对每个选中算例调用 `_read_input_cgyro_scalars(case_dir)`，返回 `[(case_name, scalars_dict), ...]`。
- `_get_flux_2d_varying_params()` → 找出所有在选中算例间 **数值不同** 的标量键。过滤条件：
  - 键名匹配 `^(Z|MASS|DENS|TEMP|DLNNDR|DLNTDR|SDLNNDR|SDLNTDR|NU)_\d+$`（物种参数）。
  - 或键名在白名单 `physical_exact` 中（核心物理/几何参数如 `KY`, `BETAE_UNIT`, `RMIN`, `Q`, `S`, `KAPPA` 等）。
  - 排除数值/求解器控制参数（如 `N_RADIAL`, `N_TIME`, `N_TOROIDAL` 等）。
- `_is_flux_2d_physical_param_key(key)` → 上述过滤规则的静态实现。
- `_refresh_flux_2d_param_choices()` → 更新 `flux_scan_xparam_combo` 的下拉值，当前选择不在列表中时自动选第一个。

**Plot Type 字符串构建**：

- `_build_effective_plot_type()` → 将 GUI 控件状态翻译成内部 `plot_type` 字符串和显示标题 `display_plot_type`。
  - Flux → `"{flux_type} Flux {xaxis_str}"`（如 `"Energy Flux vs ky"`），若勾选了分解则追加 `" (Decomp)"`。
  - Fluctuation 1D → `"Phi vs ky"` 或 `"Phi FFT"`。
  - Fluctuation 2D → `"Fluctuation 2D"` / `"Fluctuation 2D vs xt"` / `"Fluctuation 2D vs kxky"`。
  - Zonal ExB → 四种精确字符串：`"ZF ExB Shearing Rate"`、`"ZF ExB Shearing Spectrum"`、`"ZF Phi Spectrum (theta0)"`、`"ZF ExB vs gamma_lin (kx=ky)"`。
  - Energy balance → 按 mode_key 映射为 `"Energy Balance Entropy"` / `"Energy Balance ZF"` / `"Energy Balance Gamma Eff"` / `"Energy Balance Single {qty} {xaxis}"` / `"Energy Balance FULLT"` / `"Energy Balance vs 2D"`。
  - Others → `"Integration Error"` / `"Radial Correlation (rcorr_phi)"` / `"POD Parity"`。

**其他 UI 工具**：

- 类顶部 `_FLUX_*`、`_FLUC_*`、`_FFT_*`、`_ZF_*`、`_ENERGY_BALANCE_*`、`_OTHERS_*` 常量集中维护下拉菜单选项，避免字符串散落在 UI 创建和 plot-type 映射逻辑中。
- `_browse_linear_gamma_file()` → 文件选择对话框，初始目录来自 `default_share_dir()`，优先 `/data/share/{user}`。
- `_pod_theta_resolution_error_text()` / `_show_pod_theta_resolution_warning()` / `_case_theta_resolution_is_insufficient(data)` → POD 模式 theta 分辨率不足时阻止绘图并弹窗。
- `_stop_animation()` → 停止当前 `FuncAnimation`，重置控制按钮状态。
- `_enable_manual_pager(...)` / `_update_manual_pager_status()` / `prev_frame()` / `next_frame()` / `toggle_pause()` → 动画/手动翻页控制。

---

## 4. 绘图基类 Plotting（调度中枢）

### 4.1 视觉常量与样式引擎

全局常量定义在文件顶部：

```python
LINE_COLOR_PALETTE = ["#F14040", "#1A6FDF", "#37AD6B", "#B177DE", "#CC9900", ...]  # 11 色
GAMMA_LIN_COLOR = "#515151"  # 线性谱专用灰色
CONTOUR_COLOR_PALETTE = ["#08306B", "#08519C", ..., "#67000D"]  # 16 色 diverging，用于 2D contour
GRID_COLOR = "#808080"
GRID_LINEWIDTH = 0.5
GRID_ALPHA = 0.5
AXIS_BORDER_LINEWIDTH = 1.5
PLOT_LINEWIDTH = 2.0
PLOT_MARKER_SIZE = 4.0
TICK_LABEL_FONTSIZE = 14
MAIN_FONT_FONTSIZE = 14
LEGEND_FONT_FONTSIZE = 10
PREFERRED_FONT_FAMILIES = ["Arial", "DejaVu Sans", "Liberation Sans", "Noto Sans", "Helvetica"]
```

`_resolve_font_family()` 遍历系统字体，返回第一个存在的；全都不存在时回退 `"DejaVu Sans"`。

**`_apply_global_plot_color_style()`**：
- 设置 `fig.set_prop_cycle` 和 `mpl.rcParams['axes.prop_cycle']` 为 `LINE_COLOR_PALETTE`。
- 批量更新 `rcParams`：字体族、字号、标题/标签/刻度/图例字号、轴线宽、网格颜色/透明度、`axes.axisbelow=True`、`lines.linewidth=PLOT_LINEWIDTH`。

**`_apply_unified_visual_style_to_figure()`**：
- 遍历 `self.fig.axes` 中所有轴：设置 `axisbelow=True`、开启网格、调整四边 spine 线宽、刻度标签字号、所有 `ticklabel` 和 `label` 和 `title` 的字体族。
- 遍历 legend：设置文字字体。
- 遍历所有 line：统一 `linewidth=PLOT_LINEWIDTH`；如果有 marker 则设置 `markersize=PLOT_MARKER_SIZE`。
- 遍历 figure-level text（如 suptitle）设置字体。

**`_reset_figure_layout_defaults()`**：
- 关键：在每次新图前调用。`fig.clear()` **不会**重置之前 `tight_layout()` 修改的 subplot 参数，导致后续单轴图被拉伸。因此显式调用 `fig.set_layout_engine(None)` 和 `fig.subplots_adjust(left=..., right=..., ...)` 恢复 Matplotlib 默认值。

### 4.2 时间窗解析引擎

**`_get_time_indices(t)`** —— 全局唯一时间窗解析入口：

1. 将输入 `t` 转为 1D float 数组，得到 `t_min`, `t_max`。
2. `default_start = t_min + 0.5 * (t_max - t_min)`（默认后半段）。
3. 检查 `time_mode_var`：
   - `"manual start/end"` 或 `"manual range"`：读 `t_start_var` / `t_end_var`，空则 fallback 到 `default_start` / `t_max`。
   - `"full range"`：`t_start=t_min`, `t_end=t_max`。
   - `"last duration"`：读 `time_duration_var`，空则 `0.5*(t_max-t_min)`。
   - 默认 `"last percent"`：读 `time_percent_var`（默认 50），计算 `t_start = t_max - (percent/100)*(t_max-t_min)`。
4. Clamp 到 `[t_min, t_max]`。
5. 如果 `t_start >= t_end`（用户输入错误或反向），fallback 到 `default_start`。
6. `indices = np.where((t_arr >= t_start) & (t_arr <= t_end))[0]`。
7. 如果 `indices` 为空（极端稀疏网格），fallback 到后半段索引 `np.arange(start_idx, t_arr.size)`。
8. 返回 `(indices, t_start, t_end)`。

**`_compute_case_time_window(data)`**：对每个算例调用 `_get_time_indices(data.t)`。

**时间窗标注工具**：
- `_format_avg_suffix(t_start, t_end, prefix="Avg")` → `" (Avg: 10.0-20.0)"`。
- `_append_avg_suffix(label, t_start, t_end, prefix)` → 拼接。
- `_format_avg_range_from_axis(x_axis, valid_idx, prefix)` → 从轴值和索引数组反推时间范围字符串。

### 4.3 归一化与单位转换系统

**`_rho_scalar_for_norm(data, label)`**：
- 读取 `data.rho`，无效时回退 `1.0`。
- 首次使用时打印 `Info: {label} uses rho = ... for normalization.`，并通过 `_cmp_rho_norm_logged` 属性去重。

**实离子 GyroBohm 归一化（Flux 专用）**：

**`_get_flux_real_ion_norm_context(data, label)`**：
1. 检查缓存 `data._cmp_flux_real_ion_norm_ctx`，有则直接返回。
2. 解析 `input.cgyro` 标量，读取 `N_SPECIES` 以及各物种的 `Z_i`, `MASS_i`, `TEMP_i`, `DENS_i`。
3. 从 `data.z/mass/temp/dens` 数组补充缺失值。
4. **选择实离子**：`Z>0` 的物种中，`DENS` 最大的那个（`np.argmax`）。
5. 计算转换比率（基于 D-归一化 CGYRO 惯例，`m_D=1`）：
   - `vc = cs_ri/cs = sqrt(1/m_i)`
   - `rhoc = rho_s_ri/rho_s = sqrt(m_i)`
   - `gc = Gamma_GB,ri / Gamma_GB = vc * rhoc^2`
   - `qc = Q_GB,ri / Q_GB = vc * rhoc^2`（与 gc 公式相同，但物理含义不同）
6. 缓存并打印详细信息（物种号、Z、M、DENS、各比率）。
7. 失败时打印 Warning 并缓存 `{'valid': False}`。

**`_get_flux_real_ion_norm_scale(data, moment_idx, label)`**：
- 如果用户未勾选实离子归一化，返回 `1.0`。
- 如果 `moment_idx == 1`（能量通量），返回 `1.0 / qc`。
- 否则（粒子通量），返回 `1.0 / gc`。

**`_get_rhos_to_rhoe_factor(data, label)`**：
- 计算 `rho_s / rho_e` 用于 Fluctuation 2D 的空间轴转换。
- 电子质量获取优先级：
  1. `input.cgyro` 中 `MASS_{N_SPECIES}`（最后一个物种通常为电子）。
  2. `data.mass` 数组最后一个元素。
  3. 物理常数 fallback：`m_e/m_D ≈ 9.109e-28 / 3.345e-24`。
- 公式：`factor = sqrt(1.0 / m_e_norm)`（D-归一化下 `m_D=1`）。
- 缓存到 `data._cmp_rhos_to_rhoe_factor`。

**`_use_fluc2d_x_rhoe()`**：读取 `fluc2d_x_elec_var`，决定是否使用电子尺度。

### 4.4 数据加载与形状强制系统

**`_load_kxky_complex(data, label, attr_name, file_suffix, species_dependent=False)`**：
- 三层加载策略：
  1. 直接 `getattr(data, attr_name)`。
  2. 若失败且未尝试过 `getbigfield()`，调用 `_ensure_bigfield_loaded(data, label)` 后再试。
  3. 若仍失败且 `data` 有 `extract` 方法，通过 `data.extract(file_suffix, cmplx=True)` 读取二进制文件，并根据 `n_radial * theta_plot * n_n [ * n_species]` 做文件尺寸推断 reshape（Fortran order）。
- 最终通过 `_coerce_kxky_complex()` 统一形状。

**`_coerce_kxky_complex(raw, label, tag, species_dependent)`**：
- 输入可能是：复数 ndarray、`[2, ...]` 实虚部分离数组。
- 输出规范形状：
  - `species_dependent=False` → `[nr, theta, ky, t]`
  - `species_dependent=True` → `[nr, theta, species, ky, t]`
- 支持的输入变体：
  - `np.iscomplexobj(arr)` 且 ndim==4/5。
  - `arr.ndim > 0 and arr.shape[0] == 2` 的实/虚部打包格式，ndim=5/6。

**`_field_attr_suffix_from_name(field_name)`**：
- 映射：`'phi'→('kxky_phi','.cgyro.kxky_phi')`、`'apar'→('kxky_apar','.cgyro.kxky_apar')`、`'density'→('kxky_n','.cgyro.kxky_n')`、`'energy'→('kxky_e','.cgyro.kxky_e')`、`'velocity'→('kxky_v','.cgyro.kxky_v')`。

**`_load_named_kxky_complex(data, label, field_name, species_dependent=False)`**：
- 通过上述映射调用 `_load_kxky_complex`。

**`_extract_midplane_kykxt(field_complex, data, label, drop_radial0=False, species_idx=0)`**：
- 输入 `[nr, theta, ky, t]` 或 `[nr, theta, species, ky, t]`。
- 若为 5D，提取 `species_idx` 对应的物种切片。
- 使用 `_midplane_theta_index(data, n_th)` 选取 `itheta = 4*n_theta//8`（外中平面）。
- `drop_radial0=True` 时跳过第一个径向槽位（`[1:, ...]`），匹配 CGYRO 的 collect 行为。
- 输出 `[nr', ky, t]`（复数）。

**`_build_kx_axis(data, n_r, label)`**：
- 优先用 `data.kx`。常见情况：`kx.size == n_r` 直接返回；`kx.size == n_r-1` 时假设跳过了第一个槽位，返回 `kx` 和 `radial_idx=np.arange(1, n_r)`。
- 若 `kx` 缺失，用 `data.length` 计算 `dkx = 2π/length`，构建 `p_index * dkx`。
- 返回 `(kx_axis, radial_idx)`。

**`_reconstruct_x_from_kx(c_kx_t, nx=None)`**：
- 输入 `[nr, nt]` 复数数组（kx 系数随时间）。
- 通过 `np.fft.fft` 重建实空间 x 轮廓，输出 `[nt, nx]`。
- 内部使用 CGYRO 的共轭对称映射逻辑：`d[:, i] = coeff[src, :]`，其中 `i = ix if ix>=0 else ix+nx`，`src = -ix + nr//2`。
- 若 FFT 结果全零或无效，fallback 到直接 `np.matmul(coeff.T, phase)` 的慢速路径。

**`_maptoreal_fft(nr, nn, nx, ny, c)`**：
- 从 CGYRO 谱系数 `c[nr, nn]` 通过 `np.fft.irfft2` 重建实空间场。
- 关键映射：利用 `f(p,-n) = f(-p,n)*` 的厄米特性，将负 p 映射到 FFT 的正索引。
- 输出 `[nx, ny]`，乘以 `0.5` 修正半和方法。

### 4.5 1D 绘图包装

**`_plot_1d(x, y, label, plot_type)`**：
- 防御检查：`x.size==0` 或 `y.size==0` 或尺寸不匹配时直接返回。
- 仅保留 `np.isfinite(x) & np.isfinite(y)` 的点。
- 若 `plot_type` 包含 `"vs Time"`，画纯线（无 marker）；否则画 `marker='o', markersize=PLOT_MARKER_SIZE`。

### 4.6 绘图主流程

**`plot_comparison()`** —— 用户点击 **Plot** 后的总入口：

```
plot_comparison()
  ├─ _reset_plot_area()
  │   ├─ _stop_animation()
  │   ├─ 清空 energy_entropy 双轴缓存
  │   ├─ fig.clear()
  │   ├─ _reset_figure_layout_defaults()
  │   └─ _apply_global_plot_color_style()
  │   └─ ax = fig.add_subplot(111)
  ├─ _build_effective_plot_type() → (selection, plot_type, display_plot_type)
  ├─ _get_selected_case_names() → [case_name, ...]
  ├─ 若 _is_contour_like_plot(plot_type) 且 len>1：
  │   弹窗警告，仅保留第一个算例
  ├─ 若 plot_type=="POD Parity" 且 theta 分辨率不足：
  │   弹窗阻止，直接 return
  ├─ 分支分发：
  │   a) Flux + "vs 2D" → _plot_flux_vs_2d_selected_cases()
  │   b) Energy balance + "vs 2D" → _plot_energy_balance_vs_2d_selected_cases()
  │   c) plot_all_species_var 且 Flux → _plot_all_species_flux_first_case()
  │   d) 默认 → _plot_selected_cases()
  └─ _finalize_plot(selection, plot_type, display_plot_type)
```

**`_plot_selected_cases(selected_cases, plot_type)`**：
- 遍历选中算例，对每个调用 `_plot_single_case(data, case_name, plot_type)`，并捕获异常。

**`_plot_single_case(data, label, plot_type, species_override_index=None)`**：
1. `_compute_case_time_window(data)` → `(t_indices, t_start, t_end)`。
2. `_dispatch_single_case_plot(data, label, plot_type, t_indices, t_start, t_end, species_override_index)`。
3. 异常捕获并打印 traceback。

**`_dispatch_single_case_plot(...)`** —— 中央路由：
1. 先调用 `_dispatch_common_plot_families(...)` 处理模式匹配：
   - `"Frequency"` / `"Growth Rate"` → `_plot_frequency_growth()`
   - `"Flux"` in plot_type → `_plot_flux_diagnostics()`
   - `"FFT"` in plot_type → `_plot_fluctuation_fft()`
   - 返回 `True` 表示已处理，流程结束。
2. 未命中时，处理特殊歧义分支：`plot_type.startswith("Energy Balance Single ")` → `_plot_energy_balance_single_mode()`。
3. 处理通用涨落 1D：若包含 `"vs ky"` / `"vs kx"` / `"vs Time"` 且不包含 `"Flux"`，且字段名在 `["Phi","Apar","Bpar"]` 中 → `_plot_fluctuation_1d()`。
4. 最后查精确映射表 `_build_exact_plot_handler_map()`：
   - `ZF ExB Shearing Rate` → `_plot_zf_exb_shearing_rate_time_trace`
   - `ZF ExB Shearing Spectrum` → `_plot_zf_exb_shearing_rate_kx_spectrum`
   - `ZF Phi Spectrum (theta0)` → `_plot_zf_phi_kx_theta0_spectrum`
   - `ZF ExB vs gamma_lin (kx=ky)` / `ZF ExB Fig4 (kx=ky)` → `_plot_zf_exb_vs_gamma_lin_kx_equals_ky`
   - `Energy Balance Entropy` → `_plot_energy_balance_entropy`
   - `Energy Balance ZF` → `_plot_energy_balance_zf`
   - `Energy Balance Gamma Eff` → `_plot_energy_balance_gamma_eff_v3`
   - `Energy Balance Single` → `_plot_energy_balance_single_mode`
   - `Energy Balance FULLT` → `_plot_energy_balance_fullt`
   - `Fluctuation 2D` → `_plot_fluctuation_2d`
   - `Fluctuation 2D vs xt` → `_plot_xt_fluctuation_contours`
   - `Fluctuation 2D vs kxky` → `_plot_fluctuation_kxky_map_from_2d`
   - `Integration Error` → `_plot_other_error`
   - `Radial Correlation (rcorr_phi)` → `_plot_other_rcorr_phi`
   - `POD Parity` → `_plot_other_apar_pod_parity`

**`_finalize_plot(plot_type_selection, plot_type, display_plot_type)`**：
- 如果 `_energy_entropy_axes_active`（entropy vs ky 双轴模式），分别对两个轴设置 legend、grid、axisbelow，然后全局样式并 `canvas.draw()` 返回。
- 如果不是标准线图（contour-like），仅应用全局样式并返回。
- 标准线图：
  - 收集 legend handles，设置 `loc='best'` 等参数。
  - 标题：`f"Comparison: {display_plot_type}"`，若 Frequency/Growth Rate 勾选了 `/ky` 则追加 `/ky`。
  - 调用 `_apply_axis_labels(plot_type)` 根据 plot_type 自动推断 xlabel/ylabel（如果尚未设置）。
  - 应用 `log_x` / `log_y`。
  - 调用 `_apply_unified_visual_style_to_figure()` + `canvas.draw()`。

### 4.7 坐标轴标签推断

**`_apply_axis_labels(plot_type)`**：
- 如果 `self.ax.get_xlabel()` 已非空，**不覆盖**（尊重专科方法已设置的标签）。
- 按 plot_type 子串匹配推断：
  - `Frequency` → `x: $k_y \rho_s$`, `y: $\omega\ (c_s/a)$`
  - `Growth Rate` → `x: $k_y \rho_s$`, `y: $\gamma\ (c_s/a)$`
  - `ZF ExB Shearing Rate` → `x: $t\ (a/c_s)$`, `y: $\omega_{E\times B}^{ZF}\ (c_s/a)$`
  - `Flux` + `vs ky` → 根据实离子归一化决定 `$k_y \rho_{s,ri}$` 或 `$k_y \rho_s$`
  - `Flux` + `Energy` → `$Q_s/Q_{GB}$`（或 `$Q_{GB,ri}$`）
  - `Flux` + `Particle` → `$\Gamma_s/\Gamma_{GB}$`（或 `$\Gamma_{GB,ri}$`）
  - `vs Time` → `$t\ (a/c_s)$`

### 4.8 辅助诊断：Case Info 与 Input Diff

**`plot_case_info()`**：
- 无选中 → 用第一个算例；有选中 → 用第一个选中算例。
- 调用 `_plot_other_case_info(data, case_name)`，读取 `out.cgyro.info` 全文。
- 启用 **连续滚动模式**：将文本按行拆分为 `self._case_info_lines`，颜色列表 `self._case_info_line_colors`。
- `_enable_case_info_scroll()`：
  - 计算每屏可见行数 `_estimate_case_info_lines_per_view()`（基于 axes 像素高度 / 15px 每行）。
  - 绑定鼠标滚轮事件 `_on_case_info_scroll(event)`，每次滚动 3 行。
  - 启用手动翻页器，用户可用 `< Prev` / `Next >` 按钮逐屏浏览。
- `_render_case_info_window(top_line_idx)`：渲染从 `top_line_idx` 开始的可见行，标题显示 `[top+1-end/total]`。

**`plot_input_diff()`**：
- 选择规则：`>=2` 选中 → 比较所有选中；`1` 选中 → 比较选中与其他所有；`0` 选中 → 比较所有已加载。
- 总是走 `_plot_input_diff_multi_cases()`，渲染 **多算例标量差异矩阵**。
- 过程：
  1. 对每个算例解析 `input.cgyro` 为 `scalars: dict[str, float]`。
  2. 取所有共同键的交集。
  3. 筛选出数值变化范围 `>1e-12` 的键。
  4. 按行输出：`[KEY]` 下各算例的值；若某算例缺失显示 `<MISSING>`（红色）。
  5. 用 matplotlib text 渲染，不同颜色区分不同算例。

**`_plot_input_diff_two_cases(case_a, case_b)`**：
- 备选路径（历史遗留）：渲染 `difflib.unified_diff` 的彩色文本 diff（绿+/红-）。

---

## 5. 专科绘图模块详解

### 5.1 FrequencyPlotting（频率 / 增长率）

**`_plot_frequency_growth(data, label, plot_type, t_indices, t_start, t_end)`**：
- 输入：`data.freq`，形状可能是 `[2, n_ky]`（单时刻）或 `[2, n_ky, n_t]`（多时序）。
- `comp_idx = 0`（Frequency）或 `1`（Growth Rate）。
- 若 `freq.ndim == 3`：
  - 自动检测 ky 所在轴：比较 `comp.shape[0]` 和 `comp.shape[1]` 与 `n_ky`。
  - 在 `t_indices` 上取 `np.mean`；若 `t_indices` 为空，取最后时刻 `-1`。
  - 标注 `" (Avg: t0-t1)"`。
- `/ky` 归一化：若 `self.norm_ky_var.get()` 为真，用 `np.divide(y, x_safe, where=|x|>1e-12)`，零处填 `nan`。

### 5.2 FluxPlotting（通量）

**核心数据结构**：`data.ky_flux` → `[species, moment, field, ky, time]`，其中 `moment=0` 粒子，`moment=1` 能量；`field=0` Phi, `1` Apar, `2` Bpar。

**入口：`_plot_flux_diagnostics(data, label, plot_type, t_indices, t_start, t_end, species_override_index)`**：

- `moment_idx = 1 if "Energy" in plot_type else 0`。
- `is_decomp = "Decomp" in plot_type`。
- `vs_kx_estimated = "vs kx" in plot_type.lower()`。

**分支 A：vs kx (estimated)** → `_plot_flux_kx_estimated_spectrum(...)`：
- **不依赖 `ky_flux`**，直接从谱场文件反算。
- 加载 `phi = kxky_phi`, `mom_n = kxky_n`, `mom_e = kxky_e`, `mom_v = kxky_v`, `apar = kxky_apar`, `bpar = kxky_bpar`。
- 主矩：`primary = mom_e if moment_idx==1 else mom_n`。
- EM 代理：粒子模式用 `mom_v` 作为 Apar 的代理；能量模式用 `primary` 本身作为 Apar/Bpar 的代理（近似）。
- 所有数组截断到共同最小维度 `[n_r, n_th, n_ky, n_t]`。
- 通道核函数：`kernel = np.real(np.conj(moment_arr) * (1j * ky4 * field_arr * prefactor))`。
  - ES: `prefactor = -1.0`（phi 与密度的交叉相位）。
  - EM-A: `prefactor = 1.0`（Apar）。
  - EM-B: `prefactor = 1.0`（Bpar）。
- 对 `theta` 平均、对 `ky` 求和 → `flux(kx, t)`。
- 时间平均 → 按 `kx` 排序 → 只保留 `kx >= 0` 分支 → `_plot_1d()`。
- 若 `is_decomp`，叠加 ES(Phi)、EM-A(Apar)、EM-B(Bpar) 三条虚线。
- 公式面板通过 `_render_flux_kx_formula_math()` 显示推导细节。

**分支 B：vs ky_time** → 2D contour：
- `flux_sel = ky_flux[target_indices, moment_idx]` → `[species, field, ky, t]`。
- 用户指定时间窗时，对时间轴切片；否则用全部时间。
- 对 species 和 field（若未分解则仅 phi）求和 → `[ky, t]`。
- `transpose` 为 `[t, ky]` 后用 `ax.contourf(ky_plot, t_plot, z_t_ky, levels=80, cmap=...)`。
- 若 `is_decomp`：创建 `n_fields` 个纵向子图，分别画 Phi/Apar/Bpar。

**分支 C：vs ky** → 1D 线图：
- `flux_ky_t = np.sum(flux_sel, axis=(0,1))` → `[ky, t]`。
- 若 `is_decomp`：总通量实线 + 各 field 虚线（`linestyle='--', marker='x'`）。
- 若 `t_valid.size > 0`：时间平均，标注 `"Time Avg"`。
- 否则 fallback 到最后时刻。

**分支 D：vs Time** → 时间演化：
- `y_total = np.sum(flux_sel, axis=(0,1,2))` → `[t]`。
- 若 `is_decomp`：总通量 + 各 field 虚线，标注时间窗均值。
- 若 `t_valid.size > 0`：在图上画均值虚线（`linestyle='--', color=line.get_color()`）。

**分支 E：vs 2D（参数扫描）** → `_plot_flux_vs_2d_selected_cases(selected_cases, plot_type)`：
- 从 `input.cgyro` 解析各算例标量，找出变化的物理参数。
- 用户选择 X 轴参数（如 `KY`）。
- 自动分组：剩余变化参数相同的算例连成一条曲线。特殊处理：`DLNNDR_*` 系列耦合参数不会导致分组分裂。
- 对每个算例：加载 flux，时间平均，计算 `y_val`。
- 对每组数据：按 `x_val` 排序，合并重复 x（平均），用 `marker='o'` 画线。
- Y 轴标签根据能量/粒子和实离子归一化自动设置。

**`_plot_all_species_flux_first_case(selected_cases, plot_type)`**：
- 当用户勾选 `"Plot All Species (First Case Only)"` 时触发。
- 仅对第一个选中算例，遍历所有物种索引，逐个调用 `_plot_single_case(data, case_name, plot_type, species_override_index=i)`。

### 5.3 FluctuationPlotting（涨落）

**`_load_fluctuation_moment_field(data, label, moment, main_ion_policy)`**：
- 统一加载不同矩的复杂场：
  - `"Phi"/"Apar"/"Bpar"` → `_load_named_kxky_complex(..., species_dependent=False)`。
  - `"Density"/"Energy"` → `_load_named_kxky_complex(..., species_dependent=True)`，然后对选中物种求和。
  - `"Temperature"` → 需要同时加载 `Density` 和 `Energy`，按 `(2/3 * Energy - T0 * Density) / n0` 计算温度扰动。`main_ion_policy="first"` 时仅保留一个离子（温度是单物种量）。
- 输出 `[nr, theta, ky, t]` 复数数组。

**`_plot_fluctuation_1d(data, label, plot_type, t_indices, t_start, t_end)`**：
- 字段名从 `plot_type.split()[0]` 提取（Phi/Apar/Bpar）。
- 加载并 `_extract_midplane_kykxt(drop_radial0=True)` → `[nr-1, ky, t]`。
- `rho_norm = data.rho`（无效时 `1.0`），全场除以 `rho_norm`。
- 平均模式：
  - **Mean Absolute**：`y = sum(|F|)`（对 kx 或 ky 求和）。
  - **Root Mean Square**：`y = sqrt(sum(|F|^2))`。
- **vs ky**：对 kx 求和（axis 0）→ `[ky, t]` → 时间平均 → `_plot_1d(ky, y, label, plot_type)`。ylabel 为 `$\sum \langle |F/\rho_s| \rangle_t$` 或 `$\sqrt{\sum \langle |F/\rho_s|^2 \rangle_t}$`。
- **vs kx**：对 ky 求和（axis 1）→ `[kx, t]`。kx 轴优先用 `data.kx`，尺寸不匹配时动态计算 `dkx`。
- **vs Time**：分离 `n=0`（`ky_idx_0`）和 `n>0`（`ky!=0` 的掩码）两条曲线。时间窗有效时，分别计算均值并在图上画均值虚线。
- **fft**：走 `_plot_fluctuation_fft()`（见 5.4）。

**`_plot_fluctuation_2d(data, label, plot_type, t_indices, t_start, t_end)`**：
- 加载矩场 → `_extract_midplane_kykxt(drop_radial0=False)` → `[nr, ky, t]`。
- `get_frame_data(ti)`：
  - 取时间切片 `c = c_midplane[:, :, ti]`。
  - `nr = c.shape[0]`, `nn = data.n_n`, `nx = nr+1`, `ny = 2*nn-1`。
  - 调用 `_maptoreal_fft(nr, nn, nx, ny, c)` → 实空间 `[nx, ny]`。
  - 乘以 `2` 修正半和方法。
  - 构建物理坐标：`x = np.arange(nx)*2π/nx` → `xp = x/(2π) * data.length`；`y` 同理用 `ky[1]` 的周期。
  - 若 `use_x_rhoe=True`：`xp *= rho_s/rho_e`。
  - 周期延拓：在 x/y 末尾补一个周期点，满足 `fp[-1,:]=fp[0,:]`, `fp[:,-1]=fp[:,0]`。
- **动画判断**：`user_specified_time = bool(t_start_var or t_end_var)` 且 `len(t_indices) > 1` 时启动 `FuncAnimation`（interval=100ms，repeat=True）。
- 否则画静态图（取 `t_indices[-1]` 或最后时间点）。
- 启用动画控制按钮（Pause/Prev/Next）。

**`_plot_xt_fluctuation_contours(data, label, plot_type, t_indices, t_start, t_end)`**：
- `vs xt` 模式：时空等高线。
- 取 `ky=0` 的 `c_kx_t = field[:, itheta, ky_idx_0, :]`（`[nr, t]` 复数）。

**`_plot_fluctuation_kxky_map_from_2d(data, label, plot_type, t_indices, t_start, t_end)`**：
- `vs kxky` 模式：从 Fluctuation 2D 的 `Moment` 选择读取对应 `kxky_*` 场。
- 调用 `_extract_midplane_kykxt(drop_radial0=True)` 得到 `[kx, ky, t]`，对选中时间窗做 RMS 或 mean-absolute 平均。
- 最终通过 `_plot_fluctuation_kxky_map(...)` 在 `(kx, ky)` 平面上绘制二维幅值图，并登记 `X Y Z` 导出数据。
- 用户指定时间窗时切片，否则用全部时间。
- `_reconstruct_x_from_kx(c_kx_t)` → `[nt, nx]` 实空间轮廓。
- `transpose` 后 `contourf(t_plot, x_plot, f_tx, levels=80, cmap=...)`。

### 5.4 FftPlotting（FFT 谱）

**`_plot_fluctuation_fft(data, label, plot_type, t_indices)`**：
- 字段名从 `plot_type.split()[0]` 提取。
- 加载 `kxky_*` 场 → `_extract_midplane_kykxt(drop_radial0=False)` → `[nr, ky, t]`。
- 时间切片：`field_t_all = field_t_all[..., t_indices]`（若 `t_indices` 非空）。
- `dt = mean(abs(diff(t_array)))`。
- **符号约定**：默认 `freq_mult = -1.0`。尝试读取 `out.cgyro.info`，若包含 `"Ion direction: omega > 0"` 或 `"< 0"`，均保持 `-1.0`（CGYRO `e^{-iωt}` 与 FFT `e^{-i2πft}` 的匹配）。
- `omega = np.fft.fftfreq(nt, d=dt) * 2π`，沿时间轴 `np.fft.fft`。
- `np.fft.fftshift` 将零频移到中心，`omega *= freq_mult`。若 `freq_mult < 0`，再 `np.flip` 确保升序。
- 谱度量：Amplitude 模式用 `|F|` 并归一化到全局/每 ky 最大值；Power 模式用 `|F|^2`，不归一化。
- **Linear coherent** = `|Σ_kx F_kx|^2`（相干叠加）；**Nonlinear incoherent** = `Σ_kx |F_kx|^2`（非相干叠加）。

**视图分支**：
- **"Omega vs ky"**：上下两个子图。
  - 上图：`kx=0` 分量（`i_r_0 = data.n_radial // 2`）。
  - 下图：Linear 模式下相干求和；Nonlinear 模式下非相干求和。
  - 可选叠加线性频率曲线（读取 `linear_gamma_file`，画 `--` 灰色线）。叠加前检查 ky 范围重叠，无重叠时跳过并打印提示。
- **"Omega vs kx"**：上下两个子图。
  - 上图：`ky=0` 分量。
  - 下图：Linear 模式下对 ky 相干求和；Nonlinear 模式下对 ky 非相干求和。
  - 仅保留 `kx >= 0` 分支。

### 5.5 ZfPlotting（带状流）

**核心源数据提取：`_get_zf_exb_phi_kx_t(data, label)`**：
- 遵循 `collect.py:get_zfshear` 的索引约定：
  1. `ky` 索引固定为 `0`（`n=0` 的 zonal 模式）。
  2. `theta` 索引为 `itheta0 = 4 * n_theta // 8`（外中平面）。
  3. 径向切片 `[1:, ...]`（跳过第一个径向槽位）。
- 四层加载：`data.kxky_phi` → `getbigfield()` → `extract('.cgyro.kxky_phi', cmplx=True)` 文件尺寸推断。
- 形状兼容：支持 `[nr, theta, nn, nt]`（复数）、`[2, nr, theta, nn, nt]`（实/虚部分离）、`[2, nr, theta, species, nn, nt]`（含物种维）。
- `kx` 对齐：若 `kx.size == n_radial-1`（已跳过首槽），直接匹配；若 `kx.size == n_radial`，截去首元素；否则动态计算 `2π*p/length`。
- 返回 `(kx, phi_kx_t, t)`，其中 `phi_kx_t` 为 `[n_radial-1, nt]` 复数。

**`_compute_zf_exb_shearing_time_series(data, label)`**：
- `shear_rate = Σ_kx kx² * |phi_kx_t| / rho`。
- 返回 `(t, shear_rate)`。

**`_plot_zf_exb_shearing_rate_time_trace(data, label, t_indices, t_start, t_end)`**：
- 调用上述计算，画时间序列。
- 时间窗有效时：标注均值/标准差，并画均值水平虚线（`linestyle='--', color=line.get_color()`）。

**`_compute_zf_phi_theta0_kx_profile(data, label, t_indices, t_start, t_end)`**：
- 计算 `<|phi_ZF|/rho_s>_t` 作为 `y_phi`，返回 dict 含排序后的 `x(kx)`、`y_phi`、原始数组等。

**`_plot_zf_exb_shearing_rate_kx_spectrum(...)`**：
- `y = x² * y_phi`，画 1D 谱。

**`_plot_zf_phi_kx_theta0_spectrum(...)`**：
- 直接画 `y_phi` vs `kx`。

**`_plot_zf_exb_vs_gamma_lin_kx_equals_ky(...)`**：
- 读取线性谱文件（`ky, omega, gamma` 三列文本），排序并保留 `ky>=0`。
- 计算 `omega_ZF(kx) = kx² * <|phi|/rho>` 和 `ky * V_ZF_mean`，其中 `V_ZF_mean = 0.5 * sqrt(Σ|kx * <|phi|/rho>|²)`。
- 三条曲线：`γ_lin(ky)`（灰色实线+圆点）、`ω_ZF(kx)`（红色虚线+x）、`ky·V_ZF_mean`（蓝色点划线+三角）。
- **比值计算**：若用户在 `zf_gamma_lin_ky_var` 中输入了 ky 值，用 `np.interp` 在该 ky 处插值三条曲线，计算 `ω_ZF/γ_lin` 和 `ky·V_ZF/γ_lin`，并追加到 legend 标签中（如 `$[\omega_{ZF}/\gamma_{lin}=0.523]$`）。
- 在该 ky 处画垂直虚线 `axvline`（通过 `_zf_ky_marker` 属性去重，防止多算例重复画线）。
- x 轴截断到线性谱的最大 ky 范围内。

### 5.6 EnergyPlotting（能量平衡）

**核心数据**：`data.triad`，形状 `[2, species, radial, 8, n_n, n_time]`（实部/虚部分离）。
- `channel 0` = T_a（idx1）
- `channel 1` = N_a 或 T_a*（idx2）
- `channel 2` = d(T_a·δS_a)/dt（idx3）
- `channel 3` = dW_em/dt（idx4）
- `channel 4` = δS_a（idx5）
- `channel 5` = D_r,a（idx6）
- `channel 6` = D_θ,a（idx7）
- `channel 7` = D_c,a（idx8）

**`_load_triad_common(data, label, require_flux)`**：
- 若 `data.triad` 不存在，先尝试 `getbigfield()`，再尝试 `data.extract('.cgyro.triad')` 文件级回退。
- 文件回退时根据 `n_n * n_radial * n_species * 8 * 2` 推断 `n_time`，Fortran order reshape 为 `[2, species, radial, 8, n, t]`。
- 若 `require_flux=True`，还需确保 `data.ky_flux` 存在（调用 `data.getflux()`）。

**`_compute_triad_v2_terms(data, label)`**：
- 重建 cgyro_plot 的核心物理量时间序列，返回大字典：
  - `T0` = Σ_radial Re(idx1)
  - `N0` = Σ_radial Re(idx2)
  - `Ent0` = Σ_radial Re(idx3)
  - `S0` = Σ_radial Re(idx5)
  - `Wkt0` = Σ_radial Re(idx4)
  - `diss_r0/diss_th0/diss_c0` = Σ_radial Re(idx6/7/8)
  - `G_plus_Q` = 经过 `dlntdr/dlnndr/temp` 修正的驱动项
  - `Wes0` = 从 `phi(kx, ky=0, t)` 计算的静电势能 `W_es(t)`
  - `dWes_dt` = `np.gradient(Wes0, t)`
  - `dSg_dt` = `idx3_total - dWes_dt`（陀螺中心熵变率）
  - `Lz_total` = `dWes_dt - N_total`
- 时间对齐：`n_time_use = min(t.size, T0.size, Ent0.size, ..., Wes0.size)`。

**W_es 计算：`_compute_zonal_wes_time_series(data, label)`**：
- `W_es(t) = Σ_a [ n_a·z_a²/(2·T_a) · Σ_kx(|phi_kx(t)|² · (1-Γ_0,a(kx))) ]`
- `Γ_0,a = I_0(b_a)·exp(-b_a)`, `b_a = (|kx|·ρ_a)²`
- `ρ_a/ρ_ref = sqrt(|T_a·m_a|)/|z_a|`（注意这里的 `mass` 是 CGYRO 归一化质量）。

**分支 A：Entropy balance** → `_plot_energy_balance_entropy(...)`：
- 四条曲线：
  - `(T-N)^{NZ→Z}` = `T0 - N0`（黑色虚线）
  - `D_Z` = `diss_r0 + diss_th0 + diss_c0`（灰色实线）
  - `-L_{Z,total}` = `-Lz_total`（蓝色虚线）
  - `dS_Z^{(g)}/dt` = `dSg_dt`（橙色实线）
- 强制白色背景（`set_facecolor('white')`），匹配论文风格。
- 标题含 `n_e = dlnndr[1]`。

**分支 B：ZF energy balance** → `_plot_energy_balance_zf(...)`：
- 四条曲线：
  - `N^{NZ→Z}` = `N0`（黑色实线）
  - `D'_Z` = 0（灰色实线，占位）
  - `L_{Z,total}` = `Lz_total`（蓝色实线）
  - `L_{T_e}·dW_es/dt` = `dWes_dt`（红色实线）

**分支 C：Effective growth rate** → `_plot_energy_balance_gamma_eff_v3(...)`：
- 加载 triad，对物种求和 → `f[radial, 8, n, t]`。
- `idx5`（channel 4）= `ent = Re(delta S_a,k)`，时均 → `ent_avg[radial, n]`。
- `gamma_eff^NZ`：取 `idx2`（channel 1），**强制 `ky=0`（n-index 0）置零**，时均 → `-a_nz_avg / (2*ent_avg)`。
- `gamma_eff^Z`：取 `idx1 - a_nz`，时均 → `-a_z_avg / (2*ent_avg)`。
- `gamma_eff^NL = gamma_eff^NZ + gamma_eff^Z`。
- 取 radial 中心行 `row = n_radial // 2`。
- ky 轴优先用 `data.kynorm`，次之用 `data.ky`；强制 `ky>=0` 并排序。
- `_extend_to_origin(x, y)`：在 `(0,0)` 处补点，保证曲线连续性。
- **无线性文件**：画 `γ_eff^NZ`（绿方块）、`γ_eff^Z`（红三角）、`γ_eff^NL`（黑色虚线）。
- **有线性文件**：画 `γ_lin`（灰色圆点实线）、`γ_eff^NZ`、`γ_eff^Z`、`γ_eff^stable = γ_lin - γ_eff^NL`（品红虚线）。
- 若用户输入了 `n_var`（实际作为 ky 值使用），插值计算各曲线在该 ky 的值，并输出 `ratio = gamma_eff / gamma_lin` 到 legend 和 stdout。

**分支 D：Single plot** → `_plot_energy_balance_single_mode(...)`：
- 支持 quantity：`T`, `N`, `T-N`, `entropy`。
- 支持 x-axis：`vs Time`（直接时间序列）、`vs ky`（时均谱）、`vs kxky`（2D 等高线）。
- `vs ky` 且 `quantity=='entropy'` 时：走 `_plot_energy_balance_single_entropy_spectrum()`，**双轴布局**（上 ion 黑色，下 electron 红色）。
  - 使用 `self._energy_entropy_axes_active` 状态机防止重复 `fig.clear()`。
  - 每算例线型通过 `_get_case_linestyle(..., map_attr="_energy_entropy_case_style", line_styles=['-', '--', '-.', ':'])` 保持稳定。
  - 计算 `log(sum_kx <delta S_a>_t)`，对 ky 作图。
- `vs ky` 非 entropy 时：通过 `_compute_energy_balance_single_vs_ky()` 计算时均的 `T(ky)`、`N(ky)` 或 `T-N(ky)`。
- `vs kxky`：通过 `_compute_energy_balance_single_vs_kxky()` 在 `(kx, ky)` 平面上画 `contourf`。
  - 对 triad channel 0 (T) 或 channel 1 (N) 做时间平均，得到 `[radial, n_ky]` 分布。
  - radial 维度映射为 kx（优先用 `data.kx`，缺失时用 `data.length` 计算），n_n 维度映射为 ky。
  - entropy 模式不支持 vs kxky。

**分支 E：FULLT transfer map** → `_plot_energy_balance_fullt(...)`：
- 加载 `data.fullt`（若不存在则 `extract('.cgyro.fullt')` 回退）。
- 形状：`[Re/Im, kx_target, ky_target, kx_source, ky_source, t]`。
- 固定 `k_target = (0, ky_fixed)`，在 `(kx_source, ky_source)` 平面上画传递函数 `T_k^Φ(k')` 的时均等高线。

**分支 F：2D scan** → `_plot_energy_balance_vs_2d_selected_cases(selected_cases)`：
- 类似 Flux vs 2D 的扫描逻辑，但物理量为能量平衡通道。
- 计算 `N_D/S`, `T_D/S`, `N_e/S`, `T_e/S`（D=主离子，e=电子）。
- S 分母 = `sum_{a,p,n!=0} Re[idx5]`（非 zonal 的 delta S 总和）。
- 若 `n>1` 的维度缺失，fallback 到全部 n 的总和。
- 按其他变化参数分组，画四条曲线（黑方=N_D/S，黑虚方=T_D/S，红三角=N_e/S，红虚三角=T_e/S）。

### 5.7 OtherPlotting（其他诊断）

**`_plot_other_error(data, label)`**：
- 读取 `data.err1`（总误差）和 `data.err2`（RK 误差）。
- 跳过前两个点（匹配 cgyro_plot 风格）。
- 负数替换为 `nan`，强制 `set_yscale('log')`。

**`_plot_other_rcorr_phi(data, label, t_indices, t_start, t_end)`**：
- 调用 `data.kxky_select(theta_idx, field_idx, 'phi', 0)` 获取场数据。
- 对 `ky!=0` 求和得到 `y[n_r, nt]`，时间平均后 `np.roll(ave, -nx//2)`。
- `corr = np.fft.fftshift(np.fft.fft(ave))`，归一化到最大幅值。
- `delta_r = np.fft.fftshift(np.fft.fftfreq(nx)) * (2π/dk)`。
- **指数包络拟合**：`scipy.signal.hilbert` 提取包络 → `curve_fit` 拟合 `exp(-|x|/tau)` → 标注 `l_corr`。

**`_plot_other_apar_pod_parity(data, label, t_indices, t_start, t_end, field_override)`**：
- 支持 `Apar` 和 `Phi` 两种场。
- 加载 `[nr, theta, ky, t]` 复数场。
- 解析用户输入的 `kx`/`ky`：
  - `ky_sel_idx = _resolve_axis_selection_index(ky_axis, ...)`
  - `kx_sel_axis_idx = _resolve_axis_selection_index(kx_axis, ..., prefer_value=True)`（优先解释为物理值而非索引）
- **连通 kx 链**：使用全部 radial 索引（`chain_field_idx = radial_idx`）。
- **相位修复**：
  1. `np.fft.fftshift` 对 kx 维度中心化。
  2. 若 kx 单调，按 kx 排序。
  3. 调用 `_repair_connected_kx_stitching()`：
     - 计算原始平滑度评分 `raw_score = Σ|Δblock|² / Σ|block|²`。
     - 强制应用 CGYRO Eulerian 交替相位 `(-1)^p`（p 相对于 kx~0 中心），得到 `arr_alt`。
     - 计算 `alt_score`。
     - 应用 `_enforce_connected_kx_phase_continuity()`：相邻切片通过重叠内积相位对齐 `v_k *= exp(-i·arg(<v_{k-1}, v_k>))`。
     - 输出修复后的数组和描述字符串。
- **POD/SVD**：
  - 构建矩阵 `X = block.reshape(n_chain*n_theta, n_t_sel)`，去均值。
  - `np.linalg.svd(X, full_matrices=False)` → `U, S, Vh`。
  - 取前两个模式 `U[:, :2]`，reshape 为 `(n_chain, n_theta)`。
- **扩展 z/π 轴**：
  - `_build_theta_over_pi_axis(data, n_theta)`：优先用 `data.thetap`，次之用 `theta/theta_b/theta_over_pi` 等候选；若尺寸不匹配，做降采样或线性插值；最终 fallback 到 `linspace(-1,1,n_theta)`。
  - `_build_extended_z_over_pi_axis(theta_over_pi, chain_field_idx, chain_kx, p_axis)`：对每个连通 kx 块，z = theta/pi + 2*p_rel，其中 p_rel 以 kx~0 为中心。
- **奇偶分析**：
  - `_parity_even_odd_ratios(profile, z_axis)`：在对称域 `[-zsym, zsym]` 上数值积分：
    - `P_even = ∫|A(z)+A(-z)|²dz / (4∫|A(z)|²dz)`
    - `P_odd  = ∫|A(z)-A(-z)|²dz / (4∫|A(z)|²dz)`
  - 判定规则：`P_even > 0.7` → tearing；`P_odd > 0.7` → ballooning；否则 mixed。
- **能量占比**：
  - `s² = S²`（奇异值平方）。
  - `f_ball = 100*s²[ball_idx]/sum(s²)`，`f_tear` 同理，`f_mix` 为未分类模式，`f_res` 为剩余模式（>=3）。
- **布局**：1x2 子图，左=mode1，右=mode2。Real 黑实线，Imag 红虚线。x 轴范围若过大则限制在 `[-DEFAULT_POD_Z_WINDOW_PI, DEFAULT_POD_Z_WINDOW_PI]`（默认 8π）。
- **标题**：包含 ky, kx, 时间窗, 连通kx标识, stitch 方法, 能量占比。

---

## 6. 数据导出模块（CgyroDataExportMixin）

数据导出、workspace 保存/加载、linear spectrum 文件选择共用 `cgyro_comparison_bootstrap.default_share_dir()` 作为默认目录策略：优先 `/data/share/$USER`，再回退 `/data/share` 和当前工作目录。

**`transfer_bin_to_readable()`**：
- 要求用户选择输出根目录。
- 对每个选中算例：
  1. `_resolve_case_dir_for_export(data)` 定位算例目录。
  2. 遍历 `bin/` 子目录下所有 `bin.cgyro.*` 文件。
  3. 普通文件：`_extract_bin_array(data, suffix)` → `_write_flat_table(arr, out_path)`。
     - 形状未知时扁平化为 `(flat_index, value)` 或 `(flat_index, real, imag)`（复数）。
  4. `.cgyro.triad`：调用 `_export_triad_by_channel()` → 按 8 个物理通道拆分：
     - 推断形状 `[2, species, radial, 8, n, time]`（Fortran order）。
     - 对每个通道 `ch`：提取 `triad[:, :, :, ch, :, :]` → `[2, species, radial, n, time]`。
     - 调用 `_write_nd_table(arr_ch, out_path, dim_labels=["ri","species","radial","n","time"], channel_index=ch+1, channel_label=...)`。
     - 输出目录：`{case_out}/triad_channels/Sheet1_T_a.txt` ... `Sheet8_D_c_a.txt`。
     - 写入 `README_channels.txt` 说明通道映射。

**`_extract_bin_array(data, suffix)`**：
- 调用 `data.extract(suffix, cmplx=False)`，若结果非 `"null"` 且 size>0 则返回。
- 若 `suffix.startswith(".cgyro.kxky_")`，fallback 到 `cmplx=True`。

---

## 7. 事件绑定与回调链总图

```
Tkinter 主窗口
├── <<ComboboxSelected>> on plot_type_combo → update_options()
│   └── 根据 plot_type 显隐动态控件，渲染对应公式面板
├── "Plot" Button → plot_comparison()
│   ├── _reset_plot_area()
│   ├── _build_effective_plot_type()
│   ├── _get_selected_case_names()
│   ├── (分支分发)
│   └── _finalize_plot()
├── "< Prev" / "Next >" / "Pause" Buttons
│   └── prev_frame() / next_frame() / toggle_pause()
│       └── anim_update_func(frame) 或 _render_case_info_window()
├── case_listbox <Button-1> → _on_drag_start()
├── case_listbox <B1-Motion> → _on_drag_motion()
├── canvas <scroll_event> (Case Info 模式)
│   └── _on_case_info_scroll() → _render_case_info_window()
└── 菜单项
    ├── File → add_case_single / add_case_multiple / add_group
    ├── Data → transfer_bin_to_readable
    ├── Time → _set_time_last_percent / _set_time_last_duration / clear_time_range
    └── Average → update_options (触发重绘公式面板)
```

---

## 8. 扩展指南

### 8.1 新增 Plot Type（完整步骤）

假设新增 `"My Diagnostic"`：

1. **`cgyro_comparison_ui.py`**：
   - 在 `plot_types` 列表追加 `"My Diagnostic"`。
   - 在 `_init_options()` 中创建该模式需要的控件变量和实例。
   - 将新控件加入 `_hide_dynamic_options()` 的 `widgets` 列表。
   - 在 `update_options()` 中添加：`if plot_type == "My Diagnostic": self.my_combo.grid(row=row, ...); row+=1; ...`。
   - 可选：实现 `_render_my_formula_math()` 并在 `update_options()` 中调用。

2. **新增或修改专科 Mixin**（如 `cgyro_comparison_plotting_other.py`）：
   - 添加 `_plot_my_diagnostic(data, label, plot_type, t_indices, ...)` 方法。

3. **`cgyro_comparison_plotting.py`**：
   - 在 `_build_exact_plot_handler_map()` 中加入 `"My Diagnostic": lambda: self._plot_my_diagnostic(...)`。
   - 或在 `_dispatch_common_plot_families()` 中加入模式匹配。
   - 若需要自动推断坐标轴标签，在 `_apply_axis_labels()` 中添加分支。

4. **`cgyro_comparison_bootstrap.py`**（如需新增全局常量）：添加并引用。

### 8.2 修改全局样式

- 调色板：`LINE_COLOR_PALETTE` / `CONTOUR_COLOR_PALETTE`（`cgyro_comparison_plotting.py` 顶部）。
- 字体/字号：`MAIN_FONT_FONTSIZE`、`TICK_LABEL_FONTSIZE` 等常量。
- 线宽：`PLOT_LINEWIDTH`、`AXIS_BORDER_LINEWIDTH`。
- 修改后所有图自动生效（通过 `_apply_global_plot_color_style` 和 `_apply_unified_visual_style_to_figure`）。

### 8.3 新增数据导出格式

- 在 `CgyroDataExportMixin` 中新增方法。
- 调用 `_extract_bin_array()` 读取二进制数据。
- 使用 `_write_flat_table()`（1D）或 `_write_nd_table()`（N-D）输出文本。

---

## 9. 依赖与运行环境

| 依赖 | 必需 | 说明 |
|------|------|------|
| `tkinter` | 是 | GUI 框架，Python 标准库 |
| `matplotlib` | 是 | 绘图引擎，`backend_tkagg` 用于嵌入 |
| `numpy` | 是 | 所有数值计算 |
| `pygacode` | 是（运行时） | `cgyrodata` / `cgyrodata_plot` 数据读取。引导层提供 Mock 降级，但无实际数据可用 |
| `scipy` | 否 | `scipy.signal.hilbert`（rcorr 包络）、`scipy.optimize.curve_fit`（rcorr 指数拟合）。缺失时跳过拟合，不影响主流程 |

**pygacode 导入策略**：
- `cgyro_comparison_bootstrap.py` 自动向上查找 `../gacode/f2py` 和 `../gacode-master/f2py` 并加入 `sys.path`。
- 若导入失败，创建 Mock 类：`cgyrodata` 返回空数组（`ky`, `freq`, `t`, `ky_flux` 等），`getflux/getbigfield` 为 no-op。
- 这使得 GUI 可以在没有 pygacode 的环境中启动（用于 UI 调试），但无法加载真实算例。

---

## 10. 常见问题与陷阱

1. **`fig.clear()` 后布局拉伸**：`_reset_figure_layout_defaults()` 显式重置了 `subplots_adjust` 和 `layout_engine`。新增使用 `tight_layout()` 或 `subplots_adjust()` 的专科方法时，务必在 `_reset_plot_area()` 中补充恢复逻辑。

2. **动画残留**：切换 Plot Type 前必须调用 `_stop_animation()`，否则会留下定时器导致 canvas 异常。

3. **kxky 数据形状多样性**：CGYRO 输出可能以复数 ndarray、`[2, ...]` 实虚部分离、或文件级扁平数组存在。所有加载必须经过 `_coerce_kxky_complex()` 统一形状，不要假设具体维度。

4. **时间索引越界**：`_get_time_indices` 有多层 fallback，但专科方法仍需对 `t_indices` 做 `(valid_t >= 0) & (valid_t < n_t)` 过滤。

5. **物种索引不匹配**：不同算例的物种数量和顺序可能不同。`_resolve_species_indices` 通过 Z/M 精确匹配；`"Main Ion"` 依赖 `DENS>0.3` 启发式；`"All Ions"` 简单取 `Z>0`。

6. **公式面板截断**：`_draw_formula_panel` 通过像素级行高估算自动扩展 fig 高度，并通过 `_formula_font_size` 缩小长行字号。若仍截断，可增大 `base_fig_w` 或减小 `width_chars`。

7. **POD 的连通 kx 修复**：CGYRO 的 Eulerian 表示有 `(-1)^p` 交替相位。`_repair_connected_kx_stitching` 强制 unwrap 后再做 SVD，否则 POD 模式会出现虚假振荡。

8. **实离子归一化的缓存**：`_get_flux_real_ion_norm_context` 将结果缓存到 `data._cmp_flux_real_ion_norm_ctx`。若算例数据被 reload 但属性对象未变，缓存可能过期。当前实现假设 reload 后 `data` 是新对象。
