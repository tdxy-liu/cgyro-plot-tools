"""Crash-safe writes shared by workspace, settings, and data export."""

from contextlib import contextmanager
import json
import os
from pathlib import Path
import tempfile


@contextmanager
def atomic_text_file(path, *, encoding="utf-8", newline=None):
    """Replace the destination only after a complete, flushed same-directory write."""
    destination = Path(path).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding=encoding, newline=newline, delete=False,
            dir=destination.parent, prefix="." + destination.name + ".", suffix=".tmp",
        ) as handle:
            temporary = Path(handle.name)
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def write_json_atomic(path, payload, *, sort_keys=False):
    """Serialize before opening files so invalid state cannot damage a saved file."""
    content = json.dumps(payload, ensure_ascii=False, indent=2,
                         sort_keys=sort_keys, allow_nan=False) + "\n"
    with atomic_text_file(path) as handle:
        handle.write(content)
