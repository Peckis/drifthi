"""Configuration loading."""

from __future__ import annotations

import pathlib

import yaml

DEFAULT_PATH = "config.yaml"


def load_config(path: str | pathlib.Path = DEFAULT_PATH) -> dict:
    p = pathlib.Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"config file {p} not found -- run from the project directory "
            f"or pass --config /path/to/config.yaml"
        )
    with open(p, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    return cfg
