"""Fake rtl_tcp server: speaks the rtl_tcp protocol and streams uint8 noise.

Lets you test hi-observe with zero hardware:
    python tests/fake_rtl_tcp.py --port 1234 &
    hi-observe --config config.yaml --duration-h 0.01
It streams as fast as the client reads (no rate limiting), so short test runs
finish quickly.
"""

from __future__ import annotations

import argparse
import socket
import struct
import threading

import numpy as np


def serve_client(conn: socket.socket) -> None:
    conn.sendall(b"RTL0" + struct.pack(">II", 5, 29))  # pretend R820T
    rng = np.random.default_rng()
    stop = threading.Event()

    def read_cmds():
        try:
            while not stop.is_set():
                buf = conn.recv(5, socket.MSG_WAITALL)
                if len(buf) < 5:
                    break
                cmd, val = struct.unpack(">BI", buf)
                print(f"[fake-rtl] cmd 0x{cmd:02x} = {val}")
        except OSError:
            pass
        stop.set()

    threading.Thread(target=read_cmds, daemon=True).start()
    block = (rng.normal(127.5, 12.0, 1 << 18)).clip(0, 255).astype(np.uint8).tobytes()
    try:
        while not stop.is_set():
            conn.sendall(block)
    except OSError:
        pass
    finally:
        stop.set()
        conn.close()
        print("[fake-rtl] client disconnected")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=1234)
    args = ap.parse_args()
    srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    srv.bind(("127.0.0.1", args.port))
    srv.listen(1)
    print(f"[fake-rtl] listening on 127.0.0.1:{args.port}")
    while True:
        conn, addr = srv.accept()
        print(f"[fake-rtl] client {addr}")
        threading.Thread(target=serve_client, args=(conn,), daemon=True).start()


if __name__ == "__main__":
    main()
