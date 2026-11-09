"""Live GUI language switching; model tokens and stored workspaces stay stable."""

import json
import os
from pathlib import Path
import tkinter as tk
from tkinter import ttk
import weakref

try:
    from cgyro_io import write_json_atomic
except ImportError:
    from .cgyro_io import write_json_atomic


# Scientific symbols, case names, paths and exported plot data are not translated.
ZH = {
    "Plot incomplete": "部分算例未能绘制",
    "Ready": "就绪", "Plotting": "绘制中", "Loaded": "已加载", "Removed": "已移除", "Removed all": "已移除全部",
    "Open Workspace...": "打开工作区…", "Add Case Group...": "添加算例分组…", "Remove Selected": "移除选中算例",
    "Remove All...": "移除全部…", "Reload All Cases": "重新读取全部算例",
    "File": "文件", "Cases": "算例", "Data": "数据", "Plot": "绘图", "Help": "帮助",
    "Settings": "设置", "Language": "语言", "About": "关于", "Close": "关闭",
    "Cancel": "取消", "Apply": "应用", "Save": "保存", "OK": "确定", "Clear": "清空",
    "Browse": "浏览", "Browse...": "浏览…", "Select All": "全选", "Search:": "搜索：",
    "Find Next": "查找下一个", "Pause": "暂停", "Play": "播放", "‹ Prev": "‹ 上一页", "Next ›": "下一页 ›",
    "Load Workspace...": "打开工作区…", "Save Workspace...": "保存工作区…", "Exit": "退出",
    "Add Case...": "添加算例…", "Add Multiple Cases...": "批量添加算例…", "Add Group...": "添加分组…",
    "Remove Selected Cases": "移除选中算例", "Remove All Cases": "移除全部算例", "Reload Cases": "重新读取算例",
    "Lock Selection": "锁定选择", "Case Info": "算例信息", "Compare Input Files": "比较输入文件",
    "Convert Binary to Readable...": "转换二进制文件…", "Export Current Plot Data...": "导出当前图数据…",
    "Batch Export Current Plot Data...": "批量导出当前图数据…", "Plot Now": "立即绘制",
    "Clear Plot": "清空绘图", "CGYRO Transfer Explorer...": "CGYRO 转移分析…",
    "CGYRO Transfer Explorer": "CGYRO 转移分析",
    "Time Window": "时间区间", "Manual Start / End": "手动起止时间", "Manual Start/End": "手动起止时间",
    "Last Percentage...": "末段比例…", "Last Duration...": "末段时长…", "Last percent": "末段比例",
    "Last duration": "末段时长", "Full range": "完整区间", "Full Range": "完整区间", "Reset Time Window": "默认时间区间",
    "Averaging": "平均方式", "Mean Absolute": "平均绝对值", "Root Mean Square": "均方根",
    "Normalization": "归一化", "None": "无", "Max value": "最大值", "Min value": "最小值",
    "Max absolute value": "最大绝对值", "Min-Max [0, 1]": "最小—最大 [0, 1]", "Min-Max [0,1]": "最小—最大 [0,1]",
    "Custom value": "自定义值", "Custom value...": "自定义值…", "Set Custom Value...": "设置自定义值…",
    "Axis": "坐标", "Set Axis Limits...": "设置坐标范围…", "Clear Axis Limits": "默认坐标范围",
    "User Guide...": "使用说明…", "Interface Guide...": "界面与语言说明…", "Update Connection...": "更新连接…", "Check for Updates...": "检查更新…",
    "Frequency": "频率", "Growth Rate": "增长率", "Flux": "通量", "Fluctuation 1D": "一维波动",
    "Fluctuation 2D": "二维波动", "Energy balance": "能量平衡", "Zonal ExB Shearing Rate": "带状 ExB 剪切率",
    "Others": "其他诊断", "Energy": "能量", "Particle": "粒子", "Density": "密度", "Temperature": "温度",
    "v.s ky": "随 ky", "v.s kx": "随 kx", "v.s theta": "随 θ", "v.s Time": "随时间",
    "v.s ky_time": "ky–时间", "v.s kx (estimated)": "随 kx（估算）", "v.s 2D": "二维扫描",
    "vs Time": "随时间", "vs ky": "随 ky", "vs kx": "随 kx", "vs xy": "x–y 平面", "vs kxky": "kx–ky 平面",
    "Field:": "场：", "Species:": "粒子种类：", "Mode:": "模式：", "Moment:": "物理量：",
    "View:": "视图：", "X-axis:": "横轴：", "Quantity:": "物理量：", "Normalize:": "归一化：",
    "Others:": "诊断：", "Spectrum:": "谱类型：", "Time Start:": "Time start · 起点", "Time End:": "Time end · 终点",
    "Divided by ky": "除以 ky", "Decompose by Field": "按场分解", "Error bars": "误差棒",
    "Normalize to real ion (max density)": "按实离子归一化（最大密度）", "All species (first case)": "全部粒子种类（首个算例）",
    "Electron-scale spatial axes": "空间坐标按电子尺度归一化", "Log X": "X 对数", "Log Y": "Y 对数", "Log z": "Z 对数",
    "Advanced: per-case kx/ky": "逐算例选择 kx/ky", "Edit per-case kx/ky...": "设置各算例 kx/ky…",
    "FFT Settings": "FFT 设置", "Nonlinear": "非线性", "Linear": "线性", "Amplitude": "振幅", "Power": "功率",
    "Weak signals (asinh)": "弱结构增强 (asinh)", "Symmetric log": "对称对数",
    "Overlay Linear Freq": "叠加线性频率", "Linear File:": "线性数据：", "2D X Param:": "扫描参数：",
    "Entropy balance": "熵平衡", "ZF energy balance": "带状流能量平衡", "Effective growth rate": "有效增长率",
    "Single plot": "单项绘图", "FULLT transfer map": "FULLT 转移图", "2D scan": "二维扫描",
    "Normalized by max T": "按 max T 归一化", "fixed source kx:": "源 kx：", "fixed source ky:": "源 ky：",
    "fixed zonal kx:": "固定带状 kx：", "n index:": "n 索引：", "ky scan:": "ky 扫描：",
    "Theta Idx:": "θ 索引：", "ky for ratio:": "比值使用 ky：", "Case": "算例", "Min": "最小值", "Max": "最大值",
    "x lim:": "X 范围：", "y lim:": "Y 范围：", "ky (physical, blank=avg)": "ky（物理量，留空取平均）",
    "kx (physical, blank=avg)": "kx（物理量，留空取平均）",
    "Maximum positive T": "最大正 T", "Minimum negative T": "最小负 T", "Stop at repeat": "重复时停止",
    "Allow repeats": "允许重复", "Skip visited modes": "跳过已访问模式", "Smooth (classic FULLT)": "平滑（经典 FULLT）",
    "Diamonds (grid samples)": "菱形（网格采样）", "Cells": "网格单元", "Auto (1, 2, 5, 10...)": "自动（1、2、5、10…）",
    "None (original ratios)": "无（原始比例）", "Custom multipliers": "自定义倍率", "Time start": "Time start · 起点",
    "Time end": "Time end · 终点", "Step": "步骤", "Peak kx′": "极值 kx′", "Peak ky′": "极值 ky′",
    "Next kx": "下一步 kx", "Next ky": "下一步 ky", "Conjugate": "共轭",
    "Connection mode": "连接方式", "Direct connection": "直接连接", "HTTP proxy (default)": "HTTP 代理（默认）",
    "Use an existing SOCKS5 proxy": "使用现有 SOCKS5 代理", "Start an SSH dynamic SOCKS5 tunnel": "建立 SSH 动态 SOCKS5 隧道",
    "HTTP proxy URL:": "HTTP 代理地址：", "Relay host:": "中转主机：", "SSH user:": "SSH 用户：",
    "SSH port:": "SSH 端口：", "Identity file:": "密钥文件：", "Local SOCKS port:": "本地 SOCKS 端口：",
    "0 = choose automatically": "0 = 自动选择", "SSH relay": "SSH 中转",
    "Choose how the application should reach GitHub for update checks and Git pull.": "设置更新检查与 Git 拉取的连接方式。",
    "Use an SSH key or ssh-agent. Passwords are not stored by this application.": "使用 SSH 密钥或 ssh-agent；程序不保存密码。",
    "HTTP CONNECT; HTTPS certificate checks remain enabled.": "使用 HTTP CONNECT，保持 HTTPS 证书验证。",
    "Example: socks5h://127.0.0.1:1080": "示例：socks5h://127.0.0.1:1080",
    "Tip: use Ctrl/Shift to select multiple folders.": "使用 Ctrl / Shift 选择多个目录。",
    "Choose a value for each case. Leave a cell blank to average over that axis; plain numbers match the nearest physical value; use idx:n only for an explicit array index.": "为每个算例设置取值。留空对该轴取平均；数字匹配最近的物理量；仅 idx:n 表示数组索引。",
}


def get_language(root):
    root = root._root()
    if not hasattr(root, "_cgyro_language"):
        root._cgyro_language = UiLanguage(root)
    return root._cgyro_language


class UiLanguage:
    def __init__(self, root):
        self.root = root
        self.path = Path(os.environ.get("CGYRO_UI_CONFIG", str(Path.home() / ".cgyro_comparison_tool" / "ui.json")))
        self.catalog = dict(ZH)
        self.widgets = weakref.WeakSet()
        self.status, self.on_change = None, None
        try:
            saved = json.loads(self.path.read_text(encoding="utf-8")) if self.path.exists() else {}
            code = saved.get("language", "zh") if isinstance(saved, dict) else "zh"
        except (ValueError, OSError):
            code = "zh"
        self.variable = tk.StringVar(root, value=os.environ.get("CGYRO_UI_LANGUAGE", code) if code in ("zh", "en") else "en")
        self._mapped = weakref.WeakSet()
        # Only a newly opened Toplevel is translated; no global redraw/polling.
        root.bind_all("<Map>", self._mapped_window, add="+")

    @property
    def chinese(self):
        # Font detection must not silently override the user's language choice.
        return self.variable.get() == "zh"

    def pair(self, chinese, english):
        self.catalog[english] = chinese
        return chinese if self.chinese else english

    def text(self, text):
        value = str(text)
        reverse = {zh: en for en, zh in self.catalog.items()}
        english = reverse.get(value, value)
        return self.catalog.get(english, english) if self.chinese else english

    def install_menu(self, menubar, status=None, on_change=None):
        self.status, self.on_change = status, on_change
        settings = tk.Menu(menubar, tearoff=False)
        languages = tk.Menu(settings, tearoff=False)
        chinese_label = "中文" if getattr(self.root, "_cgyro_theme", {}).get("cjk") else "Chinese (zh)"
        for label, code in ((chinese_label, "zh"), ("English", "en")):
            languages.add_radiobutton(label=label, variable=self.variable, value=code, command=self.apply)
        settings.add_cascade(label="Language", menu=languages)
        help_index = next((i for i in range((menubar.index("end") or 0) + 1)
                           if menubar.type(i) not in ("separator", "tearoff") and
                           menubar.entrycget(i, "label") in ("Help", "帮助")), "end")
        menubar.insert_cascade(help_index, label="Settings", menu=settings)
        self.refresh()

    def _mapped_window(self, event):
        if isinstance(event.widget, tk.Toplevel) and event.widget not in self._mapped:
            self._mapped.add(event.widget)
            self.root.after_idle(lambda: self.translate_tree(event.widget) if event.widget.winfo_exists() else None)

    def translate_tree(self, widget):
        if isinstance(widget, (tk.Tk, tk.Toplevel)):
            title = widget.title()
            for english in ("CGYRO Transfer Explorer", "CGYRO Comparison Tool"):
                chinese = self.catalog.get(english, english)
                title = title.replace(chinese, english)
                if self.chinese:
                    title = title.replace(english, chinese)
            widget.title(title)
        if isinstance(widget, LocalizedCombobox):
            widget.retranslate()
        elif isinstance(widget, tk.Menu):
            end = widget.index("end")
            for index in range((end + 1) if end is not None else 0):
                if widget.type(index) not in ("separator", "tearoff"):
                    text = widget.entrycget(index, "label")
                    translated = self.text(text)
                    if text != translated:
                        widget.entryconfigure(index, label=translated)
        else:
            try:
                if hasattr(widget, "retranslate"):
                    widget.retranslate()
                if "text" in widget.keys() and not ("textvariable" in widget.keys() and widget.cget("textvariable")):
                    text = getattr(widget, "_form_text", str(widget.cget("text")))
                    translated = self.text(text)
                    if translated != str(widget.cget("text")):
                        widget.configure(text=translated)
                    if hasattr(widget, "_form_text"):
                        widget._form_text = translated
                if isinstance(widget, ttk.Treeview):
                    for column in widget.cget("columns"):
                        widget.heading(column, text=self.text(widget.heading(column, "text")))
                if isinstance(widget, ttk.Notebook):
                    for tab in widget.tabs():
                        widget.tab(tab, text=self.text(widget.tab(tab, "text")))
                if hasattr(widget, "title") and hasattr(widget, "_sync"):
                    widget.title = self.text(widget.title)
                    widget._sync()
                if hasattr(widget, "_cgyro_form"):
                    widget._cgyro_form._key = None
                    widget._cgyro_form.schedule()
            except tk.TclError:
                pass
        for child in widget.winfo_children():
            self.translate_tree(child)

    def refresh(self):
        self.translate_tree(self.root)
        if self.on_change:
            self.on_change()

    def apply(self):
        self.refresh()
        try:
            write_json_atomic(self.path, {"language": self.variable.get()})
            text = self.pair("语言设置已保存", "Language saved")
            if self.chinese and not getattr(self.root, "_cgyro_theme", {}).get("cjk", False):
                text = "Chinese enabled; no CJK font detected. Install Noto Sans CJK SC if characters are missing."
        except (OSError, ValueError) as exc:
            text = self.pair("本次已切换，保存失败：", "Changed for this session; could not save: ") + str(exc)
        if self.status is not None:
            self.status.set(text)


class LocalizedCombobox(ttk.Combobox):
    """Translate displayed choices without changing the model's exact tokens."""
    def __init__(self, master=None, **kwargs):
        self.language = get_language(master)
        self.model = kwargs.pop("textvariable", None)
        if self.model is None:
            self.model = tk.StringVar(master)
        self.display = tk.StringVar(master)
        self._values = tuple(kwargs.pop("values", ()))
        self._syncing = False
        super().__init__(master, textvariable=self.display, **kwargs)
        self._model_trace = self.model.trace_add("write", self._from_model)
        self._display_trace = self.display.trace_add("write", self._from_display)
        self.bind("<Destroy>", self._remove_traces, add="+")
        self.language.widgets.add(self)
        self.retranslate()

    def _remove_traces(self, event):
        if event.widget is self:
            self.model.trace_remove("write", self._model_trace)
            self.display.trace_remove("write", self._display_trace)

    def _from_model(self, *_args):
        if not self._syncing:
            self._syncing = True
            try:
                self.display.set(self.language.text(self.model.get()))
            finally:
                self._syncing = False

    def _from_display(self, *_args):
        if not self._syncing:
            self._syncing = True
            try:
                reverse = {self.language.text(value): value for value in self._values}
                self.model.set(reverse.get(self.display.get(), self.display.get()))
            finally:
                self._syncing = False

    def retranslate(self):
        self._syncing = True
        try:
            super().configure(values=[self.language.text(value) for value in self._values])
            self.display.set(self.language.text(self.model.get()))
        finally:
            self._syncing = False

    def configure(self, cnf=None, **kwargs):
        if isinstance(cnf, dict):
            kwargs = {**cnf, **kwargs}
            cnf = None
        if "values" in kwargs:
            self._values = tuple(kwargs.pop("values"))
            kwargs["values"] = [self.language.text(value) for value in self._values]
        return super().configure(cnf, **kwargs)

    config = configure

    def set(self, value):
        self.model.set(value)
