# 更新连接与 HTTP 代理（v0.2.21）

HTTP 代理默认地址：`http://47.102.120.146:18889`。
未保存连接设置的新机器默认选择 HTTP。**已有明确保存的 Direct、SOCKS5 或 SSH 配置不会被覆盖**。
此时在 `Help -> Update Connection...` 选择 **HTTP proxy (default)**，核对地址后 Save 即可。

该设置只用于 Comparison 工具的 GitHub 版本检查和 Git 更新，不代理 CGYRO 计算、文件读取或 pip，不改系统代理、全局 Git 配置或保存的 origin。

## 命令行

```sh
# 使用默认 HTTP 代理（忽略此次调用中的旧连接配置）
python cgyro_update.py --check-update --http-proxy
python cgyro_update.py --update --http-proxy

# 自定义地址；也接受不带 http:// 的 HOST:PORT
python cgyro_update.py --check-update --http-proxy http://47.102.120.146:18889

# 此次调用不使用工具内代理，沿用原来的正常网络环境
python cgyro_update.py --check-update --direct
```

这些参数也可通过 `python cgyro_comparison.py ...` 调用。
`--direct`、`--http-proxy`、`--socks5-proxy`、`--ssh-relay` 互斥。
CLI 覆盖只对当前调用生效；持久化配置请使用 GUI。

## 旧版本首次更新

旧程序尚未识别 `--http-proxy` 时，在工具的 **main 分支、无未提交改动** 的目录中执行：

```sh
git -c http.proxy=http://47.102.120.146:18889 pull --ff-only https://github.com/tdxy-liu/cgyro-plot-tools.git main
```

然后重启工具，在 Update Connection 中选择 HTTP。不要强制覆盖本机代码改动。
这个命令不会永久改变 origin；如 Git 提示认证/联网问题，按实际错误排查。

## 实现边界

- 版本查询使用 HTTP CONNECT 隧道，目标仍为 HTTPS；证书和主机名检查保持开启。
- 显式 HTTP 模式不会被 `NO_PROXY` 或其他环境代理设置静默绕过；连接失败会报错，不自动直连。
- Git 拉取保留 `--ff-only` 和干净工作区检查。GitHub SSH 地址仅在该拉取进程中通过精确 URL 重写转为 HTTPS，保存的 origin 不变。
- 不能仅设置 `http.proxy` 就认为 SSH 已经过代理。非 GitHub 的 SSH 仓库会明确提示配置为 HTTPS，或选择 SOCKS5/Direct。
- 代理 URL 支持 `http://HOST:PORT`（也支持带方括号的 IPv6）。不接受账号密码、路径、查询参数或 HTTPS-to-proxy URL；HTTP CONNECT 的目标 HTTPS 不受此限制。
- 代理服务器可以看到连接的目标地址，HTTPS 正文仍端到端加密；程序不禁用 TLS 校验。
- 原格式/压缩读取器、科学诊断与归一化均未改变。

配置文件通常位于用户目录 `.cgyro_comparison_tool/update_connection.json`，或由 `CGYRO_UPDATE_CONFIG` 指定。
本版本为 JSON 增加 `http_proxy` 字段；旧版本的 mode 字段仍保留。

## 测试

`python -m unittest discover -s tests -v` 包含离线代理回归测试。
GitHub CI 不会连接这个公共代理：自动测试覆盖配置、TLS 设置、代理不绕过、临时 URL 重写和更新安全限制。
实际代理的可达性取决于所运行机器的网络；本机实测结果不能保证所有集群均可访问该地址。
