"""Minimal pure-python client for rtl_tcp.

This deliberately avoids pyrtlsdr / librtlsdr bindings: rtl_tcp ships with
every rtl-sdr driver install (including the rtl-sdr-blog fork needed for the
RTL-SDR Blog V4) and speaks a trivial binary protocol over a socket, so
nothing here can break with Python or setuptools upgrades.
"""

from __future__ import annotations

import socket
import struct

# rtl_tcp command IDs (osmocom rtl_tcp / rtl-sdr-blog fork)
CMD_SET_FREQ = 0x01
CMD_SET_SAMPLE_RATE = 0x02
CMD_SET_GAIN_MODE = 0x03    # 0 = auto, 1 = manual
CMD_SET_GAIN = 0x04         # tenths of dB
CMD_SET_FREQ_CORRECTION = 0x05
CMD_SET_AGC_MODE = 0x08
CMD_SET_BIAS_TEE = 0x0E


def connect_autostart(host: str, port: int, spawn: bool = True,
                      timeout_s: float = 12.0):
    """Connect to rtl_tcp; if nothing listens on localhost, start it ourselves.

    Returns (RtlTcp, proc): proc is the spawned rtl_tcp process (terminate it
    when done) or None if the server was already running / remote.
    """
    import shutil
    import subprocess
    import time as _t

    proc = None
    deadline = _t.monotonic() + timeout_s
    while True:
        try:
            return RtlTcp(host, port), proc
        except OSError as err:
            local = host in ("127.0.0.1", "localhost", "::1")
            if proc is None and spawn and local:
                exe = shutil.which("rtl_tcp")
                if exe is None:
                    raise ConnectionError(
                        f"rtl_tcp is not running on {host}:{port} and the binary "
                        f"is not in PATH -- start it yourself:  "
                        f"rtl_tcp -a {host} -p {port} -b 4") from err
                print(f"[sdr] rtl_tcp not running -- starting it ({exe}) ...")
                proc = subprocess.Popen(
                    [exe, "-a", host, "-p", str(port), "-b", "4"],
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if _t.monotonic() > deadline:
                if proc is not None:
                    proc.terminate()
                raise
            _t.sleep(1.0)


class RtlTcp:
    def __init__(self, host: str = "127.0.0.1", port: int = 1234, timeout: float = 10.0):
        self.host = host
        self.port = int(port)
        self.sock = socket.create_connection((host, self.port), timeout=timeout)
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 1 << 20)
        # 12-byte header: magic "RTL0", tuner type, gain count
        hdr = self._read_exact(12)
        self.magic = hdr[:4].decode("ascii", errors="replace")
        self.tuner_type, self.tuner_gain_count = struct.unpack(">II", hdr[4:12])

    def _cmd(self, cmd_id: int, value: int) -> None:
        self.sock.sendall(struct.pack(">BI", cmd_id, value & 0xFFFFFFFF))

    def set_sample_rate(self, hz: float) -> None:
        self._cmd(CMD_SET_SAMPLE_RATE, int(round(hz)))

    def set_freq(self, hz: float) -> None:
        self._cmd(CMD_SET_FREQ, int(round(hz)))

    def set_gain_db(self, db: float) -> None:
        self._cmd(CMD_SET_GAIN_MODE, 1)
        self._cmd(CMD_SET_AGC_MODE, 0)
        self._cmd(CMD_SET_GAIN, int(round(db * 10)))

    def set_ppm(self, ppm: int) -> None:
        if ppm:
            self._cmd(CMD_SET_FREQ_CORRECTION, int(ppm))

    def set_bias_tee(self, on: bool) -> None:
        self._cmd(CMD_SET_BIAS_TEE, 1 if on else 0)

    def _read_exact(self, n: int) -> bytes:
        buf = bytearray(n)
        view = memoryview(buf)
        got = 0
        while got < n:
            r = self.sock.recv_into(view[got:], n - got)
            if r == 0:
                raise ConnectionError("rtl_tcp closed the connection")
            got += r
        return bytes(buf)

    def read_samples_raw(self, n_bytes: int) -> bytes:
        """Read n_bytes of interleaved uint8 I/Q."""
        return self._read_exact(n_bytes)

    def flush(self, n_bytes: int) -> None:
        """Discard n_bytes (e.g. right after a retune)."""
        self._read_exact(n_bytes)

    def close(self) -> None:
        try:
            self.sock.close()
        except OSError:
            pass
