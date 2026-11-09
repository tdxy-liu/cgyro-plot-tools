# CGYRO Transfer Explorer 使用说明

适用版本：v0.2.32。独立 GUI 与命令行共用同一套读取、追踪及绘图逻辑。

新增 **设置 → 语言 / Settings → Language** 即时切换中英文，并保存偏好。低频说明和设置按需展开，导出入口位于结果区下方。算例与参数可通过文件菜单保存为会话。详见 [界面说明](GUI_GUIDE.md)。

## 1. 启动与依赖

在软件文件夹中运行：

```bash
python cgyro_fullt_trace.py
```

也可在原比较工具中选择 **Plot -> CGYRO Transfer Explorer...**。它会打开独立窗口，并带入当前选中的 case；没有选中时带入全部 case。

需要 Python 3.9 或更新版本、numpy、matplotlib；GUI 另外需要 tkinter。读取 FTZ 压缩文件还需要 zstandard。可一次安装独立程序的 Python 依赖：

```bash
python -m pip install -r requirements-fullt-trace.txt
```

本程序不依赖 pygacode，命令行批量导出不需要图形桌面。tkinter 由 Python 安装程序或操作系统提供。

## 2. 数据与追踪规则

选择含原生 `bin.cgyro.fullt_asym` / `.ftz` 或 `bin.cgyro.fullt` / `.ftz` 的 case 文件夹。原格式还需要对应的 `input.cgyro`、`out.cgyro.grids` 和 `out.cgyro.time` 元数据。只读取已完整保存的时间样本。

每一步使用相同时间窗内的原生实部平均：

\[
\overline{T}(\mathbf{k},\mathbf{k}') = \frac{1}{N}\sum_{j=1}^N \operatorname{Re}T(\mathbf{k},\mathbf{k}',t_j).
\]

这是保存样本的算术平均；不取绝对值、不做 RMS，也不按时间间隔加权。若输出时间间隔不同，请根据用途判断该平均方式是否适用。某个谱点只要包含非有限样本，该点就不会参与极值选择。

1. 输入起始物理 `kx` 和固定的 `ky`，匹配最近的已保存模式；`ky` 按绝对值解释。最初匹配得到的 ky 在全部步骤中保持不变。
2. 仅在 **`ky′ = 固定 ky`** 的一行寻找极值：**Maximum positive T** 选择严格大于零的最大值；**Minimum negative T** 选择严格小于零的最小值。
3. 将选中谱点的 `kx′` 作为下一步的 kx，ky 仍使用最初匹配值，继续追踪。例如 `(0, 0.2) → (-0.17, 0.2) → (-0.34, 0.2)`。
4. **Panels / steps = 4** 表示包含起点在内的四张图，即最多三次切换后续起点。每张图的实心黑色菱形标记当前模式，所选极值的详细数值保留在坐标表和导出文件中。

全部 panel 仍显示完整的、有符号的 partner `ky′` 网格，但其他 ky 分支不参与追踪。即使 `ky′=0`、负 ky′ 或其他正 ky′ 上有更大的极值，也不会跳转过去。如果固定行中没有所需符号的有限非零 T，则停止并说明原因，不回退到其他 ky。固定 ky=0 时，全程沿 ky′=0 追踪。

**FULL_T_KX0=1 的限制**：只保存了被诊断模式 `kx=0` 的数据。如果下一步落在非零 kx，已有面板仍可查看和导出，但无法凭空继续追踪；需要计算时保存全 kx 的 FULLT 诊断。

极值链是同一时间平均诊断上的逐步选择，不是时间轨迹，也不能单凭它断言能量只沿这一条链传输。

## 3. GUI 操作

1. **Add case...** 添加一个或多个 case，可在列表中多选。
2. 在 **Trace** 页选择 `FULLT_ASYM` 或 `FULLT`，设置 **Initial kx**、**Fixed ky**、步数及正/负方向。
3. 填写 **Time start / Time end**。两项留空时取实际时间跨度的后半段。只留空 start 时，从首个时间到所填 end 的中点开始；end 留空时取最后保存时间。时间窗包括两端的样本，空窗口会报错。
4. 点击 **开始追踪 / Trace cases**。读取在后台进行，**Cancel** 可取消并保留先前结果。无需一次读入全部时间历史，同一 ky 的均值会在内存预算内复用。
5. 右侧下拉框切换 case；**多面板预览 / Stacked figure** 查看纵向排图，**追踪坐标与数值 / Trace coordinates and values** 查看每步当前模式、极值位置、T 值及下一步模式。

重复模式的处理：

| 选项 | 行为 |
| --- | --- |
| Stop at repeat（默认） | 极值回到已经访问的模式时停止，避免把循环当成新的传输链。 |
| Allow repeats | 允许循环，继续到指定面板数。 |
| Skip visited modes | 在固定 ky 行排除已访问模式，再选择剩余点中最正/最负的 T。此时可能不再是该行的全局极值，排除数会写入数据。 |

若固定 ky 行内的极值恰好并列，按较小的 `|kx′|`、再按较小的 kx 确定；并列数量写入数据。

## 4. 出图样式 / Figure style 页与共享色标

每一步一行，上下排列；右侧单一红—白—蓝色标贯穿全部 panel。共用横纵坐标范围，只在最下方标注横轴；每张图上方留白标注面板编号、带放大倍数的 T 符号和当前 kx/ky，不显示 `next kx` 或 T 数值行。当前模式使用带白边的实心黑色菱形，不显示空心圆、所选极值标记或底部标记图例。详细数值保留在坐标表和导出文件中。

色标需要科学计数时，`×10ⁿ` 显示在色标右侧上方，适用于全部色标刻度，与面板的 kx/ky 标注分开。

黑色菱形的大小随当前谱网格在图中的间距变化。缩放、还原或导出时，其中心始终对准所选模式。

| 设置 | 说明 |
| --- | --- |
| Width / Panel height | 毫米单位，默认总宽 183 mm、每个数据区高 32 mm；可选择 89 mm 单栏宽度。 |
| Font | 默认 9 pt；允许 8.5–16 pt。最终投稿时仍应按期刊和插入尺寸核对。 |
| Color scale | 默认 Weak signals (asinh)。以零为中心，过渡尺度取共用上限的 25%。例如幅值为上限 5%、50%、80% 的点分别使用单侧色带的约 9.5%、69%、90%。刻度按显示量标注。可切换 Linear 或 Symmetric log，命令行对应 `--color-scale asinh`、`linear`、`symlog`。 |
| Shared color max | 所有 panel 共用固定范围 `[-M,+M]`。留空时用**未乘面板倍数**的全部原始图的最大 `|T|` 确定范围（勾选全局归一化时为 ±1）；修改面板倍数不会重算或放大色标。手动填写可覆盖这一范围。超范围颜色饱和，双端箭头表示两端可存在越界值。 |
| Panel scaling | 默认 Auto：相对最强面板，以 `1、2、5、10、20、50…` 的倍数放大较弱面板；自动倍数不会让其超出自动色标范围。None 关闭逐面板放大，保留原始相对幅度。Custom 使用下方手动倍数。 |
| Multipliers: 1, 1, 10 | 输入或粘贴倍数时自动切换到 Custom，点击 Apply 后生效。例如 `1, 1, 10` 只将第三张乘 10，共用色标及前两张的颜色保持不变。均须为有限正数。如果提前停止追踪，多余的倍数忽略，已绘制的每一张必须有对应倍数。若想恢复自动/关闭放大，可主动选择 Auto / None。 |
| Divide every panel by ONE global max | 所有 panel 先同除一个原始全局最大 `|T|`，再应用所选面板倍数。追踪表始终报告原始 T。 |
| Map rendering | 默认 Smooth (classic FULLT)：沿用红—白—蓝配色，每个谱点中心的菱形区域保持完整颜色，再平滑渐退到单元边缘的中性色。中心区占谱点到单元边缘距离的 20%，缩短白色淡出区；全部面板使用同一规则，并采用更淡、更细的网格线。添加的节点仅是绘图几何，不作为新的物理样本，也不改变原始幅值、倍数及色标范围。Diamonds / Cells 可切回离散菱形/色块。均不对弱值做阈值删除或百分位裁剪。 |
| x / y limits | 设置所有 panel 共同的显示范围；只影响显示，极值搜索仍使用固定 ky 行上的完整 kx 范围。 |
| Context title | 默认显示 case、实际采用的起止样本时间和追踪方向，可自定义或勾选 Hide case/time title 隐藏。 |
| Colorbar label | 默认空白，可选填。 |
| Export DPI | 默认 600；支持 300–1200。 |

点击 **应用出图样式 / Apply figure style** 即可修改样式，无需重新读取数据。修改初始模式、时间窗等 Trace 参数后，需要重新点击 **开始追踪 / Trace cases** 才会重新计算。

例如第三张显示 `10 × T`，表示着色值为原始 T 的 10 倍；色标范围内的读数须除以 10 才得到原始量。如果同时勾选全局归一化，还要乘回导出元数据中的 `divisor`。放大后超出色标范围的点用端点颜色表示，精确数值请看 CSV，不能由饱和颜色反推。**有面板倍数时不能直接用颜色深浅比较原始强弱**，应结合倍数或查看原始数据。自动倍数依据每张完整谱图的有限最大绝对值计算，不随缩放坐标范围而改变；全零/全非有限值面板不放大。追踪选点始终使用原始值。状态栏显示实际应用的倍数及色标范围。

加载旧会话会保留原来显式选择的 Diamonds / Cells；若要采用新的渐变风格，在 Map rendering 选择 Smooth (classic FULLT)。

PDF 和 SVG 中的文字、坐标轴、注释可编辑；高密度数据层按选定 DPI 栅格化，以控制文件大小。多 case 分别生成图组，各图组内部共用色标；若希望不同 case 之间也用相同范围，请填写相同的 **Shared color max**。

## 5. 导出与复现

选择输出父目录，然后点击 **Export figure + trace data...**，会导出当前所有成功追踪的 case。默认包含 PDF、SVG、PNG；可勾选同时导出 TIFF。同名目录已经存在时增加 `_2`、`_3` 等后缀，不覆盖原结果。

每个图组包含：

- `trace.pdf` / `trace.svg` / `trace.png` / 可选 `trace.tiff`：纵向排列的图片。
- `trace.csv`：每一步的当前模式、选中 partner、原始平均 T、下一步模式、并列数及排除数。为兼容旧文件格式保留的 `conjugated` 列在固定 ky 追踪中为 false。
- `plot_data.csv`：全部原生网格点的原始 `T_mean`、着色使用的 `display_value` 和该面板的 `display_multiplier`。满足 `display_value = T_mean / divisor × display_multiplier`，不导出绘图用的边缘零值节点。
- `maps.npz`：所有面板的原始均值矩阵，矩阵维度为 `[ky′, kx′]`，同时保存坐标和采用的时间点。
- `trace.json`：版本、读取布局、输入文件签名、时间索引、参数、实际固定的 ky、共享色标及其参考方式、各面板越界饱和点数量、实际倍数、渲染方式及停止原因。
- `README.txt`：数据解释。

**File -> Save session... / Load session...** 保存和恢复 case 路径及选项；会话文件不复制模拟数据，也不保存预览图，加载后重新 Trace 即可。

## 6. 命令行批量出图

三张图、追踪最正的 T：

```bash
python cgyro_fullt_trace.py --case /path/to/case --kx 0 --ky 0.2 --steps 3 --direction positive --start 300.5 --end 600 --output ./FULLT_exports
```

四张图、追踪最负的 T：

```bash
python cgyro_fullt_trace.py --case /path/to/case --kx 0 --ky 0.2 --steps 4 --direction negative --start 300.5 --end 600 --output ./FULLT_exports
```

多个 case，共用指定色标范围，输出单栏宽度与 TIFF：

```bash
python cgyro_fullt_trace.py --case /path/to/case_A --case /path/to/case_B --ky 0.2 --steps 4 --color-max 0.02 --width-mm 89 --formats pdf svg png tiff --output ./FULLT_exports
```

只将第三张放大 10 倍：

```bash
python cgyro_fullt_trace.py --case /path/to/case --ky 0.2 --steps 3 --direction negative --panel-multipliers 1 1 10 --output ./FULLT_exports
```

`--panel-multipliers` 自动选择 Custom 模式。默认 `--panel-scaling auto` 自动放大较弱面板；`--panel-scaling none` 保留原始比例。可追加 `--diagnostic fullt`、`--cycles allow`、`--cycles skip`、`--normalize`、`--glyph smooth`、`--glyph cells`、`--xlim -0.6 0.6`、`--ylim -0.3 0.6` 等选项。`--gui` 可带相同参数打开预设好的 GUI。查看完整参数：

```bash
python cgyro_fullt_trace.py --help
```

Windows 路径含空格时请用双引号包围。没有 DISPLAY 的计算节点使用 `--case ... --output ...` 即可；无需启动 GUI。

## 7. 常见情况

- **不足三/四张图**：查看停止原因；常见原因是重复模式、固定 ky 行没有同号的非零 T，或只有 kx=0 的诊断。只有允许重复才会继续循环。
- **No native diagnostic**：检查所选 Diagnostic 与文件名是否一致，以及是否已将原生 FULLT 文件复制到 case 目录。
- **压缩文件缺少解码器**：在运行程序的同一个 Python 环境安装 `requirements-compression.txt`；程序不会在导出期间联网安装依赖。
- **原格式文件在读取期间被改写**：本次结果会拒绝导出以避免混合不同快照，待写入完成再运行。FTZ 使用已提交的快照，后续追加数据不会改变本次选取的时间点。
- **图中其他 ky 分支的颜色更强，却没有追踪它**：极值选择限定在最初匹配的固定 ky 行；其他分支仍显示，但不参与径向追踪。
