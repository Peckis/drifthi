"""HI4PI survey access: tile download + declination-strip cache.

HI4PI (HI4PI Collaboration 2016, A&A 594, A116) is distributed by CDS as
20x20 deg FITS cubes. In the EQ2000/CAR set the tiles are named
CAR_<row><col>.fits with rows A..I centered at dec -80, -60, ... +80 and
columns 01..18 centered at RA 10, 30, ... 350 deg. Each tile is 266x266
pixels of 5 arcmin and 933 VRAD channels of 1.288 km/s (+-600 km/s), in
brightness-temperature Kelvin. (Verified against the tile headers.)

For an 80 cm dish (beam ~17 deg FWHM) full resolution is pointless, so this
module downloads only the tiles overlapping the observed declination strip
and reduces them into a single small cache file on a 0.5 deg / 2.58 km/s
grid: strip_cache.npz {ra, dec, v, cube[ra, dec, v]}.
"""

from __future__ import annotations

import argparse
import pathlib
import urllib.request

import numpy as np

from .config import load_config

ROWS = "ABCDEFGHI"          # dec centers -80 .. +80
TILE_HALF_DEG = 11.1        # 266 px * 5' / 2, includes the ~1 deg overlap


def row_center_dec(row: str) -> float:
    return -80.0 + 20.0 * ROWS.index(row)


def _circ_dist(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


def tiles_for_dec_range(dec_min: float, dec_max: float,
                        ra_min: float | None = None,
                        ra_max: float | None = None) -> list[str]:
    """Tiles overlapping a dec band; optionally restricted to an RA interval
    (interval may wrap through 360)."""
    if ra_min is not None and ra_max is not None:
        width = (ra_max - ra_min) % 360.0 or 360.0
        ra_mid = (ra_min + width / 2.0) % 360.0
        half = width / 2.0
    names = []
    for r in ROWS:
        c = row_center_dec(r)
        if c + TILE_HALF_DEG < dec_min or c - TILE_HALF_DEG > dec_max:
            continue
        for c2 in range(1, 19):
            col_center = 10.0 + 20.0 * (c2 - 1)
            if ra_min is not None and ra_max is not None and \
                    _circ_dist(col_center, ra_mid) > half + TILE_HALF_DEG:
                continue
            names.append(f"CAR_{r}{c2:02d}.fits")
    return names


def download_tile(name: str, base_url: str, dest_dir: pathlib.Path) -> pathlib.Path:
    dest = dest_dir / name
    if dest.exists() and dest.stat().st_size > 1 << 20:
        return dest
    url = base_url.rstrip("/") + "/" + name
    tmp = dest.with_suffix(".part")
    print(f"[hi4pi] downloading {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "drifthi/0.1"})
    with urllib.request.urlopen(req, timeout=120) as resp, open(tmp, "wb") as fh:
        total = int(resp.headers.get("Content-Length") or 0)
        got = 0
        while True:
            block = resp.read(1 << 20)
            if not block:
                break
            fh.write(block)
            got += len(block)
            if total:
                print(f"\r[hi4pi]   {got/1e6:7.1f} / {total/1e6:.1f} MB", end="")
        print()
    tmp.replace(dest)
    return dest


def _tile_axes(hdr):
    """Per-pixel RA, dec, velocity arrays from a tile header (CAR, CRVAL2=0)."""
    i = np.arange(hdr["NAXIS1"])
    j = np.arange(hdr["NAXIS2"])
    k = np.arange(hdr["NAXIS3"])
    ra = (hdr["CRVAL1"] + (i + 1 - hdr["CRPIX1"]) * hdr["CDELT1"]) % 360.0
    dec = (j + 1 - hdr["CRPIX2"]) * hdr["CDELT2"]
    v_kms = (k + 1 - hdr["CRPIX3"]) * hdr["CDELT3"] / 1000.0  # VRAD m/s -> km/s
    return ra, dec, v_kms


def build_strip_cache(cfg: dict, dec_min: float, dec_max: float,
                      out_path: pathlib.Path | None = None,
                      keep_tiles: bool = True,
                      ra_min: float | None = None,
                      ra_max: float | None = None) -> pathlib.Path:
    from astropy.io import fits

    h = cfg["hi4pi"]
    cache_dir = pathlib.Path(h["cache_dir"])
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_path or cache_dir / "strip_cache.npz"

    sky_bin = float(h["sky_bin_deg"])
    v_bin_ch = int(h["v_bin_ch"])
    v_max = float(h["v_max_kms"])

    ra_grid = np.arange(0.0, 360.0, sky_bin) + sky_bin / 2
    dec_grid = np.arange(dec_min, dec_max, sky_bin) + sky_bin / 2
    nra, ndec = len(ra_grid), len(dec_grid)

    tiles = tiles_for_dec_range(dec_min, dec_max, ra_min, ra_max)
    print(f"[hi4pi] dec strip [{dec_min:.1f}, {dec_max:.1f}] deg -> {len(tiles)} tiles "
          f"(~{0.252*len(tiles):.1f} GB download, cached in {cache_dir})")

    sum_cube = cnt_cube = v_out = None
    for n, name in enumerate(tiles, 1):
        path = download_tile(name, h["base_url"], cache_dir)
        print(f"[hi4pi] [{n}/{len(tiles)}] binning {name}")
        with fits.open(path, memmap=True) as hdul:
            hdr = hdul[0].header
            ra_px, dec_px, v_px = _tile_axes(hdr)
            ksel = np.where(np.abs(v_px) <= v_max)[0]
            k0, k1 = ksel[0], ksel[-1] + 1
            k1 = k0 + ((k1 - k0) // v_bin_ch) * v_bin_ch
            if v_out is None:
                v_out = v_px[k0:k1].reshape(-1, v_bin_ch).mean(axis=1)
                sum_cube = np.zeros((nra, ndec, len(v_out)), dtype=np.float64)
                cnt_cube = np.zeros((nra, ndec), dtype=np.int64)

            jsel = np.where((dec_px >= dec_min) & (dec_px < dec_max))[0]
            if len(jsel) == 0:
                continue
            data = hdul[0].data[k0:k1, jsel[0]:jsel[-1] + 1, :]  # (nv, nj, ni)
            data = np.nan_to_num(np.asarray(data, dtype=np.float32), nan=0.0)
            nv = (k1 - k0) // v_bin_ch
            data = data.reshape(nv, v_bin_ch, data.shape[1], data.shape[2]).mean(axis=1)

            ri = np.clip((ra_px / sky_bin).astype(int), 0, nra - 1)
            dj = np.clip(((dec_px[jsel] - dec_min) / sky_bin).astype(int), 0, ndec - 1)
            flat_idx = (ri[None, :] * ndec + dj[:, None]).ravel()  # (nj*ni,)
            nbins = nra * ndec
            cnt_cube += np.bincount(flat_idx, minlength=nbins).reshape(nra, ndec)
            for kk in range(nv):
                sum_cube[:, :, kk] += np.bincount(
                    flat_idx, weights=data[kk].ravel().astype(np.float64),
                    minlength=nbins).reshape(nra, ndec)
        if not keep_tiles:
            path.unlink()

    with np.errstate(invalid="ignore", divide="ignore"):
        cube = (sum_cube / cnt_cube[:, :, None]).astype(np.float32)
    cube[cnt_cube == 0] = np.nan
    np.savez_compressed(out_path, ra=ra_grid.astype(np.float32),
                        dec=dec_grid.astype(np.float32),
                        v=v_out.astype(np.float32), cube=cube)
    print(f"[hi4pi] strip cache written: {out_path} "
          f"({out_path.stat().st_size/1e6:.1f} MB, {nra}x{ndec}x{len(v_out)})")
    return out_path


def load_strip(path: str | pathlib.Path):
    z = np.load(path)
    return z["ra"], z["dec"], z["v"], z["cube"]


def beam_smooth_strip(ra, dec, cube, beam_fwhm_deg: float):
    """Smooth the strip cube to the telescope beam (Gaussian approx).

    RA axis wraps; the RA kernel width is scaled by 1/cos(dec_center) so the
    kernel is round on the sky at the strip center. Adequate for a wide beam
    on a +-15 deg strip. NaNs handled by normalized convolution.
    """
    from scipy.ndimage import gaussian_filter

    sky_bin = float(ra[1] - ra[0])
    dec0 = float(np.mean(dec))
    sigma_deg = beam_fwhm_deg / 2.3548
    sig_ra = sigma_deg / sky_bin / max(np.cos(np.radians(dec0)), 0.2)
    sig_dec = sigma_deg / sky_bin

    w = np.isfinite(cube).astype(np.float32)
    filled = np.nan_to_num(cube, nan=0.0)
    # weight by cos(dec) so the average is a true solid-angle average
    cosd = np.cos(np.radians(dec)).astype(np.float32)[None, :, None]
    num = gaussian_filter(filled * w * cosd, sigma=(sig_ra, sig_dec, 0),
                          mode=("wrap", "nearest", "nearest"))
    den = gaussian_filter(w * cosd, sigma=(sig_ra, sig_dec, 0),
                          mode=("wrap", "nearest", "nearest"))
    with np.errstate(invalid="ignore", divide="ignore"):
        sm = num / den
    return sm.astype(np.float32)


def make_track_interpolator(ra_grid, dec_grid, v_grid, smoothed_cube):
    """RA-wrap-safe interpolator over the beam-smoothed strip cube."""
    from scipy.interpolate import RegularGridInterpolator

    pad = 4
    ra_p = np.concatenate([ra_grid[-pad:] - 360.0, ra_grid, ra_grid[:pad] + 360.0])
    cube_p = np.concatenate([smoothed_cube[-pad:], smoothed_cube,
                             smoothed_cube[:pad]], axis=0)
    itp = RegularGridInterpolator((ra_p, dec_grid, v_grid), cube_p,
                                  bounds_error=False, fill_value=np.nan)

    def eval_track(ra_t, dec_t, v_out) -> np.ndarray:
        """Model spectra (n_t, n_v_out) along a (ra_t, dec_t) track."""
        ra_t = np.mod(np.asarray(ra_t, dtype=float), 360.0)
        dec_t = np.asarray(dec_t, dtype=float)
        nt, nv = len(ra_t), len(v_out)
        pts = np.empty((nt, nv, 3))
        pts[:, :, 0] = ra_t[:, None]
        pts[:, :, 1] = dec_t[:, None]  # out-of-strip decs yield NaN spectra
        pts[:, :, 2] = v_out[None, :]
        return itp(pts.reshape(-1, 3)).reshape(nt, nv)

    return eval_track


def track_spectra(ra_grid, dec_grid, v_grid, smoothed_cube,
                  ra_t, dec_t, v_out) -> np.ndarray:
    """One-shot convenience wrapper around make_track_interpolator."""
    return make_track_interpolator(ra_grid, dec_grid, v_grid,
                                   smoothed_cube)(ra_t, dec_t, v_out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Download HI4PI tiles for the observed dec strip and build the cache")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--dec-min", type=float, default=None)
    ap.add_argument("--dec-max", type=float, default=None)
    ap.add_argument("--margin-deg", type=float, default=None,
                    help="half-width of the strip around the pointing dec "
                         "(default: beam FWHM + search range)")
    ap.add_argument("--delete-tiles", action="store_true",
                    help="delete each 252 MB tile after binning it")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    if args.dec_min is None or args.dec_max is None:
        import time
        from .velocity import pointing_radec
        _, dec0 = pointing_radec(np.array([time.time()]), cfg)
        margin = args.margin_deg
        if margin is None:
            margin = cfg["pointing"]["beam_fwhm_deg"] + cfg["compare"]["search_deg"]
        dec_min = float(dec0[0] - margin)
        dec_max = float(dec0[0] + margin)
        print(f"[hi4pi] pointing dec ~ {dec0[0]:+.1f} deg, strip +- {margin:.0f} deg")
    else:
        dec_min, dec_max = args.dec_min, args.dec_max
    dec_min = max(dec_min, -90.0)
    dec_max = min(dec_max, 90.0)
    build_strip_cache(cfg, dec_min, dec_max, keep_tiles=not args.delete_tiles)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
