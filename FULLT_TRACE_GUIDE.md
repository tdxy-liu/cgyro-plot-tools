# FULLT Trace Studio 使用说明

适用版本：v0.2.26。独立 GUI 与命令行共用同一套读取、追踪及绘图逻辑。

## 1. 启动与依赖

在软件文件夹中运行：

```bash
python cgyro_fullt_trace.py
```

也可在原比较工具中选择 **Plot -> FULLT Trace Studio...**。它会打开独立窗口，并带入当前选中的 case；没有选中时带入全部 case。

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

1. 输入起始物理 `kx`、`ky`，匹配最近的已保存模式；`ky` 按绝对值解释。
2. **Maximum positive T** 选择严格大于零的最大值；**Minimum negative T** 选择严格小于零的最小值。
3. 找到该值的 partner 坐标 `(kx′, ky′)`，将其作为下一步模式，重复此过程。
4. **Panels / steps = 4** 表示包含起点在内的四张图，即最多三次切换后续起点。最后一张仍会标记其极值。

全部 panel 保留有符号的 partner `ky′`。若选中 `ky′ < 0`，下一步采用存储半平面的对应模式 **`(-kx′, -ky′)`**，同时反转两个分量。图头 `(*)`、坐标表和导出数据都会记录此变换；它不是单独把 `ky′` 取绝对值。若对应模式没有保存，程序停止，不将它强行匹配到别的模式。

**FULL_T_KX0=1 的限制**：只保存了被诊断模式 `kx=0` 的数据。如果下一步落在非零 kx，已有面板仍可查看和导出，但无法凭空继续追踪；需要计算时保存全 kx 的 FULLT 诊断。

极值链是同一时间平均诊断上的逐步选择，不是时间轨迹，也不能单凭它断言能量只沿这一条链传输。

## 3. GUI 操作

1. **Add case...** 添加一个或多个 case，可在列表中多选。
2. 在 **Trace** 页选择 `FULLT_ASYM` 或 `FULLT`，设置起始 kx、ky、步数及正/负方向。
3. 填写 **Time start / Time end**。两项留空时取实际时间跨度的后半段。只留空 start 时，从首个时间到所填 end 的中点开始；end 留空时取最后保存时间。时间窗包括两端的样本，空窗口会报错。
4. 点击 **Trace selected cases**。读取在后台进行，**Cancel** 可取消并保留先前结果。无需一次读入全部时间历史，同一 ky 的均值会在内存预算内复用。
5. 右侧下拉框切换 case；**Stacked figure** 查看纵向排图，**Trace coordinates** 查看每步当前模式、极值位置、T 值及下一步模式。

重复模式的处理：

| 选项 | 行为 |
| --- | --- |
| Stop at repeat（默认） | 极值回到已经访问的模式时停止，避免把循环当成新的传输链。 |
| Allow repeats | 允许循环，继续到指定面板数。 |
| Skip visited modes | 排除已访问模式，再选择剩余点中最正/最负的 T。此时可能不再是原图的全局极值，排除数会写入数据。 |

若极值恰好并列，优先选择无需半平面变换的点，再按较小的 `|ky′|`、`|kx′|` 和 kx 顺序确定；并列数量写入数据。

## 4. Publication 页与共享色标

每一步一行，上下排列；右侧单一色标贯穿全部 panel。共用横纵坐标范围，只在最下方标注横轴；每张图上方留白显示当前 kx/ky、下一步 kx 和选中的原始 T 值。空心圆为当前模式，空心菱形为所选极值。

| 设置 | 说明 |
| --- | --- |
| Width / Panel height | 毫米单位，默认总宽 183 mm、每个数据区高 32 mm；可选择 89 mm 单栏宽度。 |
| Font | 默认 9 pt；允许 8.5–16 pt。最终投稿时仍应按期刊和插入尺寸核对。 |
| Shared color max | 默认自动使用全部 panel 的最大 `|T|`，色标范围为 `[-M,+M]`，中心为零。手动填写更小范围时显示越界色标箭头。 |
| Divide every panel by ONE global max | 所有 panel 同除一个全局最大 `|T|`，保留相对强弱；不对每张图单独放大。图头和追踪表仍报告原始 T。 |
| Map marks | 默认 Diamonds：每个菱形对应一个原生谱点，颜色为该点均值。Cells 使用网格色块。两者都没有插值平滑或伪造零值节点。 |
| x / y limits | 设置所有 panel 共同的显示范围；只影响显示，极值搜索仍使用完整网格。 |
| Context title | 默认显示 case、实际采用的起止样本时间和追踪方向，可自定义或勾选 Hide case/time title 隐藏。 |
| Colorbar label | 默认空白，可选填。 |
| Export DPI | 默认 600；支持 300–1200。 |

点击 **Apply publication style to preview** 即可修改样式，无需重新读取数据。修改初始模式、时间窗等 Trace 参数后，需要重新点击 **Trace selected cases** 才会重新计算。

PDF 和 SVG 中的文字、坐标轴、注释可编辑；高密度数据层按选定 DPI 栅格化，以控制文件大小。多 case 分别生成图组，各图组内部共用色标；若希望不同 case 之间也用相同范围，请填写相同的 **Shared color max**。

## 5. 导出与复现

选择输出父目录，然后点击 **Export figure + trace data...**，会导出当前所有成功追踪的 case。默认包含 PDF、SVG、PNG；可勾选同时导出 TIFF。同名目录已经存在时增加 `_2`、`_3` 等后缀，不覆盖原结果。

每个图组包含：

- `trace.pdf` / `trace.svg` / `trace.png` / 可选 `trace.tiff`：纵向排列的图片。
- `trace.csv`：每一步的当前模式、选中 partner、原始平均 T、下一步模式、是否变换半平面、并列数及排除数。
- `plot_data.csv`：全部网格点的原始 `T_mean` 与实际着色使用的 `display_value`。
- `maps.npz`：所有面板的原始均值矩阵，矩阵维度为 `[ky′, kx′]`，同时保存坐标和采用的时间点。
- `trace.json`：版本、读取布局、输入文件签名、时间索引、参数、共享色标及停止原因。
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

可追加 `--diagnostic fullt`、`--cycles allow`、`--cycles skip`、`--normalize`、`--glyph cells`、`--xlim -0.6 0.6`、`--ylim -0.3 0.6` 等选项。`--gui` 可带相同参数打开预设好的 GUI。查看完整参数：

```bash
python cgyro_fullt_trace.py --help
```

Windows 路径含空格时请用双引号包围。没有 DISPLAY 的计算节点使用 `--case ... --output ...` 即可；无需启动 GUI。

## 7. 常见情况

- **不足三/四张图**：查看停止原因；常见原因是重复模式、没有同号的非零 T、只有 kx=0 的诊断，或下一步落在未保存的 Nyquist 边界。只有允许重复才会继续循环。
- **No native diagnostic**：检查所选 Diagnostic 与文件名是否一致，以及是否已将原生 FULLT 文件复制到 case 目录。
- **压缩文件缺少解码器**：在运行程序的同一个 Python 环境安装 `requirements-compression.txt`；程序不会在导出期间联网安装依赖。
- **原格式文件在读取期间被改写**：本次结果会拒绝导出以避免混合不同快照，待写入完成再运行。FTZ 使用已提交的快照，后续追加数据不会改变本次选取的时间点。
- **负 ky′ 上的极值与下一张的 kx 符号不同**：这是同时翻转 `(kx′,ky′)` 到已保存半平面的结果，详见图头 `(*)` 和坐标表。
