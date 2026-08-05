# Update connection configuration

The update checker can reach GitHub directly, through an existing SOCKS5
proxy, or through an SSH relay that the application starts temporarily.

## GUI

Open `Help -> Update Connection...` and choose:

- `Direct connection` for the normal network path.
- `Use an existing SOCKS5 proxy` when `ssh -D` has already been started.
- `Start an SSH dynamic SOCKS5 tunnel` to let the application run `ssh -D`
  for each update check or Git pull.

The SSH mode requires an SSH key or an ssh-agent. The application never saves
an SSH password. A local SOCKS port of `0` selects a free local port.

## Command line

Start a tunnel automatically for a version check or update:

```text
python cgyro_update.py --check-update --ssh-relay relay.example.cn --ssh-user your_user --ssh-key ~/.ssh/id_ed25519
python cgyro_update.py --update --ssh-relay relay.example.cn --ssh-user your_user --ssh-key ~/.ssh/id_ed25519
```

If a SOCKS5 tunnel is already running:

```text
python cgyro_update.py --check-update --socks5-proxy socks5h://127.0.0.1:1080
python cgyro_update.py --update --socks5-proxy socks5h://127.0.0.1:1080
```

After saving the settings in the GUI, no proxy arguments are needed. The
settings are stored per user in:

```text
~/.cgyro_comparison_tool/update_connection.json
```

Use `--direct` to ignore the saved settings for one command.
