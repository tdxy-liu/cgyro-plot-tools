# CGYRO Comparison Tool

这是一个基于 Python 的图形用户界面 (GUI) 工具，用于比较和可视化 CGYRO 模拟结果。它允许用户加载多个模拟案例，并绘制各种物理量随 $k_y$ 或时间的演化图，以及 2D 涨落动画。

## 功能特性

### 最新更新
- **File -> Save workspace / Load workspace**: 可将当前加载的案例列表、选中案例和主要绘图选项保存为 JSON workspace，并在之后恢复同一绘图环境。默认保存/加载目录优先使用 `/data/share/$USER`。
- **Data -> Save current plot data**: 可将当前图中的曲线或二维图数据导出为 Origin 友好的竖排文本表。单数据集格式为 `X  Y  Z`，多数据集格式为 `X  Y  Z  Dataset`。
- **统一导出目录**: readable bin export、plot data export 和 workspace 文件默认从 `/data/share/$USER` 开始选择，若该目录不存在则回退到 `/data/share` 或当前工作目录。
- **时间窗模式**: 支持 Manual Start/End、Full range、Last percent 和 Last duration，便于快速切换全时段、后半段统计或固定长度平均窗口。
- **一致的案例配色**: Fluctuation/Phi 相关一维图中，同一个 case 的 `n=0` 与 `n>0` 使用一致颜色族，并通过线型区分，方便多案例对比。
- **Energy balance / FULLT 可视化增强**: Energy balance 增加 `vs kxky` 映射，FULLT transfer map 使用固定 source ky，显示目标 `(kx, ky)` 平面上的传输分布，并支持实部、虚部和绝对值选择。

### 1. 数据加载与管理
- **Add Case**: 选择单个 CGYRO 模拟目录加载。
- **Add Group**: 选择父目录，批量加载其下的所有有效 CGYRO 模拟子目录。
- **Remove/Remove All**: 从列表中移除已加载的案例。
- **Species Detection**: 自动检测加载案例中包含的粒子种类（如 Electron, Deuterium, Tritium, Impurities 等），并在下拉菜单中列出。
- **Workspace Save/Load**: 保存或恢复当前案例列表、选中状态和绘图控制参数。

### 2. 绘图类型 (Plot Types)

工具支持多种绘图模式，通过下拉菜单选择：

*   **Frequency**: 模的实频 $\omega$ vs $k_y$。
*   **Growth Rate**: 模的增长率 $\gamma$ vs $k_y$。
*   **Energy Flux vs ky**: 能量通量谱 $Q/Q_{GB}$ vs $k_y$。
*   **Energy Flux vs Time**: 能量通量随时间演化 $Q(t)$。
*   **Particle Flux vs ky**: 粒子通量谱 $\Gamma/\Gamma_{GB}$ vs $k_y$。
*   **Particle Flux vs Time**: 粒子通量随时间演化 $\Gamma(t)$。
*   **Phi FFT**: 电势 $\phi$ 的频谱分析（$\omega$ vs $k_y$ 等高线图）。
*   **Fluctuation 2D**: 2D 湍流涨落动画（$x$ vs $y$）。
*   **Energy Balance**: 支持 zonal/non-zonal energy balance、单量曲线、`vs kxky` 映射和 FULLT transfer map。

### 3. 选项与控制

*   **Species**: 选择要分析的粒子种类。
*   **Divide by ky**: 将 Frequency 或 Growth Rate 除以 $k_y$ 进行归一化显示。
*   **Plot All Species (First Case Only)**:
    *   仅适用于 Flux 相关绘图。
    *   勾选后，将忽略 "Species" 选择，直接绘制第一个选中案例中**所有**粒子种类的通量曲线，方便进行多粒子对比。
*   **Time Range (Start/End)**:
    *   指定绘图的时间窗口。
    *   对于 "vs ky" 图，取该时间段内的平均值。
    *   对于 "vs Time" 图，在该时间段内计算平均值和标准差。
    *   对于 "Fluctuation 2D"，播放该时间段内的动画。
*   **Moment**:
    *   仅用于 "Fluctuation 2D" 模式。
    *   可选择绘制的物理量：`Phi` (电势), `Density` (密度), `Energy` (能量), `Temperature` (温度), `Apar`, `Bpar`。
*   **Data Export**:
    *   `Save current plot data` 将当前图像数据导出为 tab 分隔的 `X Y Z` 文本，适合直接导入 Origin。
    *   对多条曲线或多个二维数据集，导出文件会额外包含 `Dataset` 列用于区分来源。

## 使用说明

1.  **启动程序**: 运行 `python cgyro_comparison.py`。
2.  **加载数据**: 点击 "Add Case" 或 "Add Group" 加载数据目录。
3.  **选择绘图**:
    *   在左侧列表中选择一个或多个案例。
    *   在 "Plot Type" 下拉菜单中选择所需的图表类型。
4.  **配置参数**:
    *   根据需要设置时间范围、粒子种类等选项。
    *   若要查看所有粒子的通量对比，勾选 "Plot All Species"。
5.  **绘图**: 点击 "Plot" 按钮生成图像。
6.  **动画**: 若选择了 "Fluctuation 2D" 并指定了时间范围，程序将自动播放动画。

## 注意事项

*   **Fluctuation 2D** 和 **Phi FFT** 模式一次只能处理一个案例。如果选择了多个案例，程序将默认只处理第一个。
*   **Temperature** 涨落是根据密度 ($\delta n$) 和能量 ($\delta e$) 涨落计算得出的：$\delta T = \frac{2}{3}\delta e - T_0 \delta n$。
*   程序会自动尝试识别粒子名称（如 Electron, Deuterium），如果无法识别，将显示为 `Ion (Z=..., M=...)`。

## 依赖库

*   `tkinter`
*   `matplotlib`
*   `numpy`
*   `pygacode` (GACODE 的 Python 接口)
