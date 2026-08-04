# CGYRO Comparison Tool 用户说明书

适用版本：0.2.11 及以上  
说明书文件：`USER_GUIDE.md`

本工具用于加载、比较和可视化多个 CGYRO 模拟案例。说明书面向日常使用者，重点介绍界面操作、绘图流程、工作区保存、数据导出、更新和常见问题。

在程序中可通过 `Help -> User Guide...` 或按 `F1` 随时打开本说明书。

## 1. 启动前准备

### 1.1 运行环境

程序需要以下组件：

- Python 3
- tkinter
- numpy
- matplotlib
- pygacode（GACODE 的 Python 接口）

在 OMFIT 或已配置好 GACODE Python 环境的终端中进入程序目录：

```text
cd /path/to/cgyro_comparison_tool
python cgyro_comparison.py
```

### 1.2 案例目录

需要选择能够被 pygacode 正常读取的 CGYRO 输出目录。若案例缺少必要输出文件、计算未完成或文件权限不足，加载可能失败。

## 2. 五分钟快速上手

1. 点击左侧 `Add Case` 加载一个案例，或使用 `Add Multiple` / `Add Group` 批量加载。
2. 在 Cases 列表中选择需要比较的一个或多个 case。
3. 在 `Plot type` 中选择绘图类别。
4. 在中间 `Options` 区设置时间、物种、坐标轴和归一化选项。
5. 点击左下方固定的 `Plot` 按钮。
6. 使用右侧 Matplotlib 工具栏缩放、平移或保存图片。
7. 需要下次继续时，使用 `File -> Save Workspace...` 保存当前工作状态。

## 3. 界面布局

### 3.1 顶部菜单

- `File`：打开、保存 workspace，退出程序。
- `Cases`：添加、删除、重载 case。
- `Data`：转换 CGYRO binary 文件，导出当前绘图数据。
- `Plot`：立即绘图、清图、设置时间窗、平均方式和坐标轴范围。
- `Help`：打开本说明书、检查更新、查看版本信息。

### 3.2 左侧栏

左侧栏分为三个区域：

- 顶部固定 Cases：case 列表、选择锁定和 case 管理按钮始终可见。
- 中间 Plot setup：仅此区域滚动；不同 Plot type 会显示不同选项。
- 底部固定 Actions：Plot、Case info、Diff input、Clear plot 和动画/翻页控制始终可见。

case 列表拥有独立滚动条。鼠标在列表上滚动时，不会带动 Plot setup 区域一起滚动。

### 3.3 右侧绘图区

右侧为 Matplotlib 绘图区。底部工具栏通常包含：

- Home：恢复初始视图。
- Back / Forward：切换历史视图。
- Pan：平移。
- Zoom：框选缩放。
- Configure subplots：调整子图边距。
- Save：保存当前图像。

窗口最下方状态栏显示 case 数量、选择数量、当前绘图类型和操作状态。

## 4. Case 管理与选择

### 4.1 添加案例

- `Add Case`：直接选择一个 CGYRO case 目录。
- `Add Multiple`：先选择父目录，再从检测到的有效子目录中多选。
- `Add Group`：扫描父目录并加载其中所有有效 case 子目录。

同名案例会自动添加后缀，避免覆盖已加载案例。

### 4.2 删除、重载和排序

- `Remove`：删除 Cases 列表中选中的案例。
- `Remove All`：确认后移除所有案例。
- `Reload`：从原目录重新读取已加载案例。
- 拖拽 Cases 列表中的名称可以调整顺序；绘图和图例通常按该顺序生成。

### 4.3 Case 选择规则

- 单击选择一个 case。
- 使用 Ctrl 或 Shift 选择多个 case。
- 默认情况下，如果没有显式选择 case，绘图会使用全部已加载 case。

### 4.4 Lock selection

勾选 `Lock selection` 后，程序会记住当前显式选择：

- 切换 Plot type、下拉框或输入框时不会退回“全部 case”。
- 后续重新选择 case 时，锁定内容会同步更新。
- 删除 case 后，失效的锁定项会自动清理。
- 锁定状态和选择会随 workspace 保存、恢复。

如果希望恢复“未选择即使用全部 case”，关闭 Lock selection 并清空列表选择。

## 5. 公共绘图选项

### 5.1 时间窗

`Time Start` 和 `Time End` 用于设置统计或显示时间范围。还可以使用 `Plot -> Time Window`：

- `Manual Start / End`：使用左侧输入框。
- `Last Percentage...`：使用数据末尾指定百分比。
- `Last Duration...`：使用末尾指定时间长度。
- `Full Range`：使用完整时间范围。
- `Reset Time Window`：清除时间限制。

时间窗为空或无有效采样点时，部分绘图会回退到最后一个时间切片；具体定义可查看 Options 中的公式说明。

### 5.2 坐标与平均

- `Log X` / `Log Y`：使用对数坐标。
- `Plot -> Averaging`：为涨落量选择 Mean Absolute 或 Root Mean Square。
- `Plot -> Axis -> Set Axis Limits...`：设置 kx、ky 等坐标范围。
- `Clear Axis Limits`：恢复自动范围。

## 6. Plot type 说明

### 6.1 Frequency

绘制实频率随 ky 的变化。`Divided by ky` 可显示除以 ky 后的结果。适合多案例线性结果比较。

### 6.2 Growth Rate

绘制增长率随 ky 的变化。选项与 Frequency 类似，也支持 `Divided by ky`。

### 6.3 Flux

主要选项包括：

- Flux type：Energy 或 Particle。
- X-axis：`v.s ky`、`v.s kx (estimated)`、`v.s Time`、`v.s ky_time`、`v.s 2D`。
- Species：选择物种。
- Decompose by Field：按不同场贡献拆分。
- Plot All Species：仅对第一个 case 绘制全部物种。
- Normalized by real ion：按真实主离子参数归一化。

默认 GyroBohm 标签使用氘参考标记，例如 `Q_GB,D`。

### 6.4 Fluctuation 1D

可选择字段 `Phi`、`Apar`、`Bpar`，以及：

- `v.s ky`
- `v.s kx`
- `v.s Time`
- `v.s theta`
- `fft`

对 `v.s ky`、`v.s kx` 和 `v.s theta`：

- 固定方向输入留空表示对该方向取平均。
- `idx:n` 表示按数组索引选择，例如 `idx:0`。
- `val:x` 表示选择最接近物理值 x 的模式，例如 `val:0.24`。
- `Normalize by max value` 将最终有限数据按最大值归一化。

#### Advanced：每个 case 单独设置 kx / ky

1. 选择 Fluctuation 1D。
2. 选择 `v.s ky`、`v.s kx` 或 `v.s theta`。
3. 勾选 `Advanced: per-case kx/ky`。
4. 点击编辑按钮。
5. 为每个 case 分别输入 kx 和 ky；留空表示该方向取平均。
6. 点击 Apply。

Advanced 配置会随 workspace 保存和恢复。

#### FFT

FFT 模式支持 Linear / Nonlinear、Omega vs ky / Omega vs kx，以及 Amplitude / Power。线性对比文件可通过文件选择按钮加载。

### 6.5 Fluctuation 2D

- `vs xy`：实空间动画。
- `vs xt`：时空图。
- `vs kxky`：kx-ky 频谱平面。

Moment 可选择 Phi、Density、Energy、Temperature、Apar、Bpar。Density、Energy、Temperature 通常还需要选择 Species。

动画生成后，使用底部 `< Prev`、`Pause`、`Next >` 控制。非动画的多页文本结果也会复用这些按钮翻页。

### 6.6 Energy balance

包含以下模式：

- Entropy balance
- ZF energy balance
- Effective growth rate
- Single plot
- FULLT transfer map
- 2D scan

不同模式会显示对应的 species、n/ky、quantity、x-axis、normalization 和 source kx/ky 选项。Options 下方的公式面板给出当前模式使用的计算定义。

### 6.7 Zonal ExB Shearing Rate

支持：

- vs Time
- vs kx
- phi vs kx(theta=0)
- vs gamma_lin

`vs gamma_lin` 需要选择线性 gamma 文件，并可指定用于比值的 ky。

### 6.8 Others

当前包含 Error、rcorr_phi 和 POD_parity。POD_parity 需要足够的 theta 分辨率；分辨率不足时程序会给出警告。

## 7. Workspace 保存与恢复

### 7.1 保存

使用 `File -> Save Workspace...` 或 Ctrl+S。Workspace 是 JSON 文件，包含：

- case 路径和顺序
- 当前选择和 Lock selection
- 主要 Plot options
- Advanced 每 case kx/ky
- 窗口大小和左右分栏位置
- 当前是否需要恢复后重绘

### 7.2 打开

使用 `File -> Open Workspace...` 或 Ctrl+O。程序会重新读取 case 路径并恢复界面状态。

Workspace 保存的是 case 路径而不是完整模拟数据。如果 case 被移动、删除或当前机器无法访问该路径，相应案例无法恢复。

## 8. 数据查看与导出

### 8.1 Case info

点击 `Case info` 查看所选 case 的 CGYRO 信息。内容较长时可使用底部翻页控制。

### 8.2 Diff input

点击 `Diff input` 比较多个 case 的输入差异。建议先显式选择需要比较的 case。

### 8.3 导出当前绘图数据

使用 `Data -> Export Current Plot Data...` 或 Ctrl+Shift+E。

- 一维曲线导出为适合 Origin 的文本列。
- 二维数据通常导出 X、Y、Z。
- 多数据集会增加 Dataset 列区分来源。

必须先成功绘图，才有可导出的当前图数据。

### 8.4 Binary 转 readable

使用 `Data -> Convert Binary to Readable...`，将选定 case 的 `bin.cgyro.*` 数据转换为可读文本。默认输出目录优先为 `/data/share/$USER`，不可用时回退到其他可写目录。

## 9. 在线与命令行更新

### 9.1 GUI 更新

使用 `Help -> Check for Updates...`。发现新版本后，程序会询问是否执行 Git fast-forward 更新。更新成功后程序自动重启，并尝试恢复更新前的 workspace 和图形。

### 9.2 命令行检查

```text
python cgyro_comparison.py --check-update
```

也可以直接运行：

```text
python cgyro_update.py --check-update
```

### 9.3 命令行更新

```text
python cgyro_comparison.py --update
```

或：

```text
python cgyro_update.py --update
```

命令行更新要求程序目录是 Git checkout、当前分支与目标分支一致，并且更新可以安全 fast-forward。

## 10. 快捷键

- F1：打开用户说明书。
- Ctrl+O：打开 workspace。
- Ctrl+S：保存 workspace。
- Ctrl+Q：退出。
- Ctrl+N：添加单个 case。
- Ctrl+Shift+N：添加多个 case。
- Delete：Cases 列表有焦点时删除选中 case。
- F5：重载 case。
- Ctrl+Shift+E：导出当前绘图数据。
- Ctrl+Enter：立即绘图。

## 11. 常见问题

### 11.1 切换选项后又绘制了全部 case

先显式选择 case，再勾选 `Lock selection`。标题下方会显示已选择的 case 数量。

### 11.2 滚动 Cases 时 Plot setup 也滚动

当前版本中两个滚动区域相互独立。若仍出现该问题，请确认正在运行的版本不低于 0.2.10。

### 11.3 git pull 提示 local changes would be overwritten

说明本地文件有未提交修改。先备份或提交本地改动；也可以在确认后使用 Git stash。不要直接删除不确定来源的修改。

### 11.4 SSH 连接 GitHub 22 端口超时

可将 origin 改为 HTTPS，或配置可用的 SSH/SOCKS 中转。更改前先执行 `git remote -v` 确认当前远端。

### 11.5 cgyro_update.py 出现 SyntaxError

确认使用 Python 3：

```text
python --version
python3 cgyro_update.py --check-update
```

### 11.6 找不到 pygacode

请在正确的 GACODE/OMFIT 环境中启动，并检查 pygacode 的 Python 路径。普通系统 Python 未必包含该模块。

### 11.7 绘图为空或报数据维度错误

- 确认 case 计算完成且输出文件完整。
- 点击 Reload 重新读取。
- 缩小到单个 case 排查。
- 清空时间范围后重试。
- 检查所选 Moment、Species、kx/ky 是否存在于该 case。

## 12. 说明书与开发文档

- 本文件 `USER_GUIDE.md`：日常操作说明，Help 菜单直接读取。
- `README.md`：功能概览和版本更新摘要。
- `DEVELOPMENT_GUIDE.md`：开发与维护说明。
- `cgyro_comparison_user_manual.tex`：原有 LaTeX 操作手册源码。
- `cgyro_comparison_manual.tex`：原有 LaTeX 模块说明源码。

如果 Help 窗口提示找不到说明书，请确认 `USER_GUIDE.md` 与 `cgyro_comparison_ui.py` 位于同一目录。
