"""Raw-session data format: a directory of .npz chunks plus meta.json.

Each chunk holds N frequency-switching cycles:
    t_mid   (N,)        unix time at cycle midpoint
    on      (N, nfft)   mean |FFT|^2 power spectrum at freq_on (fftshifted)
    off     (N, nfft)   mean |FFT|^2 power spectrum at freq_off (fftshifted)
    n_on    (N,)        number of FFTs averaged into `on`
    n_off   (N,)        number of FFTs averaged into `off`
"""

from __future__ import annotations

import json
import pathlib
import time

import numpy as np


def new_session_dir(root: str | pathlib.Path, tag: str = "obs") -> pathlib.Path:
    stamp = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    d = pathlib.Path(root) / f"{stamp}_{tag}"
    d.mkdir(parents=True, exist_ok=False)
    return d


def write_meta(session_dir: str | pathlib.Path, meta: dict) -> None:
    with open(pathlib.Path(session_dir) / "meta.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)


def read_meta(session_dir: str | pathlib.Path) -> dict:
    with open(pathlib.Path(session_dir) / "meta.json", "r", encoding="utf-8") as fh:
        return json.load(fh)


def write_chunk(session_dir, idx: int, t_mid, on, off, n_on, n_off) -> pathlib.Path:
    p = pathlib.Path(session_dir) / f"chunk_{idx:06d}.npz"
    tmp = p.with_suffix(".tmp.npz")
    np.savez_compressed(
        tmp,
        t_mid=np.asarray(t_mid, dtype=np.float64),
        on=np.asarray(on, dtype=np.float32),
        off=np.asarray(off, dtype=np.float32),
        n_on=np.asarray(n_on, dtype=np.int64),
        n_off=np.asarray(n_off, dtype=np.int64),
    )
    tmp.replace(p)  # atomic-ish: never leave a half-written chunk_*.npz
    return p


def load_session(session_dir: str | pathlib.Path):
    """Return (meta, t_mid, on, off) with cycles concatenated and time-sorted."""
    d = pathlib.Path(session_dir)
    meta = read_meta(d)
    chunks = sorted(d.glob("chunk_*.npz"))
    if not chunks:
        raise FileNotFoundError(f"no chunk_*.npz files in {d}")
    ts, ons, offs = [], [], []
    for c in chunks:
        with np.load(c) as z:
            ts.append(z["t_mid"])
            ons.append(z["on"])
            offs.append(z["off"])
    t = np.concatenate(ts)
    on = np.concatenate(ons, axis=0)
    off = np.concatenate(offs, axis=0)
    order = np.argsort(t)
    return meta, t[order], on[order], off[order]
