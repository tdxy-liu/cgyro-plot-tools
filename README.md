# cgyro-plot-tools
# CGYRO Comparison Tool (CGYRO 比较工具)


这是一个基于 Python 的图形用户界面 (GUI) 工具，旨在方便用户加载、比较和可视化 CGYRO 模拟结果。它构建于 `tkinter` 和 `matplotlib` 之上，并利用 `pygacode` 库来读取 CGYRO 数据。


## 主要功能


*   **多案例管理**: 可以加载多个 CGYRO 模拟案例目录，并在列表中进行管理。

*   **拖放排序**: 支持通过拖放操作重新排列已加载的案例，从而调整绘图顺序。

*   **灵活的绘图选项**: 提供多种物理量的可视化，包括频率、增长率、通量、电势等。

*   **动态选项界面**: 根据选择的绘图类型，自动显示相关的配置选项（如粒子种类、归一化选项、FFT 设置等）。

*   **动画支持**: 对于随时间变化的二维数据（如 Fluctuation 2D），支持播放动画以观察演化过程。

*   **数据导出**: 虽然主要用于可视化，但通过 matplotlib 的工具栏可以保存生成的图表。


## 依赖项


运行此工具需要以下 Python 库：


*   `tkinter` (通常随 Python 安装)

*   `matplotlib`

*   `numpy`

*   `pygacode` (GACODE 套件的一部分，用于读取 CGYRO 数据)


## 使用方法


### 启动工具


在终端中运行以下命令启动工具：


```bash

python cgyro_comparison.py

```


### 操作流程


1.  **加载数据**:

    *   点击 "Add Case" 选择单个 CGYRO 运行目录（包含 `input.cgyro` 或 `out.cgyro.freq` 等文件）。

    *   点击 "Add Group" 选择一个包含多个子目录的父目录，批量加载所有子案例。


2.  **选择绘图类型 (Plot Type)**:

    在左侧面板的下拉菜单中选择主要绘图类型：

    *   **Frequency**: 线性频率。

    *   **Growth Rate**: 线性增长率。

    *   **Flux**: 能量或粒子通量。

    *   **Fluctuation 1D**: 场的一维涨落。

    *   **Fluctuation 2D**: 二维涨落云图。


3.  **配置选项 (Options)**:

    根据选择的绘图类型，配置下方的选项：

    *   **Time Start / Time End**: 设置时间平均或动画的时间范围（永久置顶显示）。

    *   **Frequency / Growth Rate**: 可勾选 "Divided by ky" 进行归一化。

    *   **Flux**: 选择 "Energy" 或 "Particle" 通量，以及 "v.s ky" 或 "v.s Time"。选择特定粒子种类或 "Main Ion (D+T)"。

    *   **Fluctuation 1D**: 选择 "v.s ky"、"v.s Time" 或 "fft"。若选择 "fft"，可进一步配置 FFT 模式（线性/非线性）和视图（Omega vs ky/kx）。

    *   **Fluctuation 2D**: 选择物理量（Phi, Density, Energy 等）。对于特定物理量，需要选择对应的粒子种类。


4.  **绘图**:

    *   在案例列表中选择一个或多个案例（按住 Ctrl 或 Shift 进行多选）。注意：某些绘图类型（如 FFT 和 2D 涨落）仅支持单案例绘图。

    *   点击 "Plot" 按钮生成图表。

    *   点击 "Clear Plot" 清除当前绘图。


5.  **动画控制**:

    *   如果绘制了 Fluctuation 2D 并且指定了时间范围，可以使用下方的 "Prev", "Pause/Play", "Next" 按钮控制动画。


## 绘图类型详解


*   **Frequency / Growth Rate**: 绘制线性模拟的实频和增长率随 $k_y$ 的变化。支持除以 $k_y$ 的归一化显示。

*   **Flux**:

    *   **Energy/Particle Flux**: 绘制能量或粒子通量。

    *   **vs ky**: 显示通量谱。

    *   **vs Time**: 显示通量随时间的演化（包含均值和标准差统计）。

*   **Phi**:

    *   **vs ky**: 电势涨落幅值随 $k_y$ 的分布。

    *   **vs Time**: 电势涨落幅值随时间的演化。

    *   **Phi FFT**: 对电势进行时域 FFT 分析，显示频率谱 ($\omega$ vs $k_y$ 或 $\omega$ vs $k_x$)。

*   **Fluctuation 2D**: 在垂直于磁力线的平面 ($x, y$) 上绘制物理量的空间涨落结构。支持 Phi, Density, Energy, Temperature, Apar, Bpar 等物理量。


## 注意事项


*   请确保 `pygacode` 已正确安装并在 `PYTHONPATH` 中，或者脚本能够通过相对路径找到它。

*   FFT 分析依赖于足够的时间点数据，请确保选择的时间窗口内有足够的数据点。

*   部分高级绘图功能（如 Fluctuation 2D）需要读取较大的输出文件，加载可能会有短暂延迟。

