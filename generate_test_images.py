#!/usr/bin/env python3
"""Generate PALETTE test images (BMP + PPM). Requires only numpy.

Usage: python make_palette_tests.py [output_dir]
Output: <output_dir>/bmp/*.bmp, <output_dir>/ppm/*.ppm, <output_dir>/bmp_edge/*.bmp
"""

import colorsys
import os
import struct
import sys

import numpy as np

S = 400  # side length of the square test images


def hex_rgb(h):
    h = h.lstrip("#")
    return (
        np.array([int(h[i : i + 2], 16) for i in (0, 2, 4)], dtype=np.float64) / 255.0
    )


def to_u8(a):
    return (np.clip(a, 0.0, 1.0) * 255).round().astype(np.uint8)


# ---------------------------------------------------------------- writers


def write_ppm(path, arr):
    """Binary PPM (P6). arr is an (h, w, 3) uint8 RGB array."""
    h, w, _ = arr.shape
    with open(path, "wb") as f:
        f.write(f"P6\n{w} {h}\n255\n".encode("ascii"))
        f.write(arr.tobytes())


def write_bmp(path, arr, bits=24, top_down=False):
    """Uncompressed BMP. arr is (h, w, 3) uint8 RGB. bits is 24 or 32."""
    h, w, _ = arr.shape
    bpp = bits // 8
    stride = ((w * bpp + 3) // 4) * 4  # rows padded to a multiple of 4
    rows = range(h) if top_down else range(h - 1, -1, -1)
    pix = bytearray()
    for y in rows:
        bgr = arr[y][:, ::-1]  # RGB -> BGR
        if bits == 32:
            pad = np.zeros((w, 1), dtype=np.uint8)
            bgr = np.concatenate([bgr, pad], axis=1)
        row = bgr.tobytes()
        pix += row + b"\x00" * (stride - len(row))
    offset = 14 + 40
    file_hdr = struct.pack("<2sIHHI", b"BM", offset + len(pix), 0, 0, offset)
    info_hdr = struct.pack(
        "<IiiHHIIiiII",
        40,
        w,
        -h if top_down else h,
        1,
        bits,
        0,
        len(pix),
        2835,
        2835,
        0,
        0,
    )
    with open(path, "wb") as f:
        f.write(file_hdr + info_hdr + pix)


# ---------------------------------------------------------------- images


def four_flat_colors():
    h = S // 2
    a = np.zeros((S, S, 3))
    a[:h, :h] = hex_rgb("#d7263d")
    a[:h, h:] = hex_rgb("#1b998b")
    a[h:, :h] = hex_rgb("#2e294e")
    a[h:, h:] = hex_rgb("#f6ae2d")
    return to_u8(a)


def hue_sweep():
    a = np.zeros((S, S, 3))
    for x in range(S):
        a[:, x] = colorsys.hsv_to_rgb(x / S, 0.85, 0.9)
    return to_u8(a)


def lightness_ramp():
    a = np.zeros((S, S, 3))
    for y in range(S):
        a[y, :] = hex_rgb("#3b82f6") * (0.15 + 0.85 * y / (S - 1))
    return to_u8(a)


def near_identical_pairs():
    a = np.zeros((S, S, 3))
    a[:] = hex_rgb("#8a8a8a")
    a[50:180, 40:170] = hex_rgb("#c0392b")
    a[50:180, 230:360] = hex_rgb("#c4402f")  # almost the same red
    a[230:360, 40:170] = hex_rgb("#2874a6")
    a[230:360, 230:360] = hex_rgb("#2a78a8")  # almost the same blue
    return to_u8(a)


def rare_colors():
    a = np.zeros((S, S, 3))
    a[:] = hex_rgb("#f3ead8")
    a[25:41, 25:41] = hex_rgb("#ff0066")  # 256 px
    a[340:354, 330:344] = hex_rgb("#00c853")  # 196 px
    a[190:198, 190:198] = hex_rgb("#6200ea")  # 64 px
    return to_u8(a)


def noisy_five_blobs(seed=42):
    rng = np.random.default_rng(seed)
    centers = [
        hex_rgb(c) for c in ("#e74c3c", "#f1c40f", "#2ecc71", "#3498db", "#9b59b6")
    ]
    labels = rng.integers(0, 5, size=(S, S))
    a = np.zeros((S, S, 3))
    for i, c in enumerate(centers):
        a[labels == i] = c
    a += rng.normal(0, 0.04, a.shape)
    return to_u8(a)


def odd_size_four_colors(w=401, h=301):
    """Odd width so 24-bit BMP rows need padding (401 * 3 = 1203 bytes)."""
    a = np.zeros((h, w, 3), dtype=np.uint8)
    a[: h // 2, : w // 2] = (215, 38, 61)
    a[: h // 2, w // 2 :] = (27, 153, 139)
    a[h // 2 :, : w // 2] = (46, 41, 78)
    a[h // 2 :, w // 2 :] = (246, 174, 45)
    return a


# ---------------------------------------------------------------- main


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "palette_tests"
    dirs = {k: os.path.join(root, k) for k in ("bmp", "ppm", "bmp_edge")}
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    main_images = {
        "1_four_flat_colors": four_flat_colors,
        "2_hue_sweep": hue_sweep,
        "3_lightness_ramp": lightness_ramp,
        "4_near_identical_pairs": near_identical_pairs,
        "5_rare_colors": rare_colors,
        "6_noisy_five_blobs": noisy_five_blobs,
    }
    for name, fn in main_images.items():
        img = fn()
        write_bmp(os.path.join(dirs["bmp"], name + ".bmp"), img)
        write_ppm(os.path.join(dirs["ppm"], name + ".ppm"), img)

    odd = odd_size_four_colors()
    write_bmp(os.path.join(dirs["bmp_edge"], "odd_width_24bit_bottom_up.bmp"), odd)
    write_bmp(
        os.path.join(dirs["bmp_edge"], "odd_width_24bit_top_down.bmp"),
        odd,
        top_down=True,
    )
    write_bmp(os.path.join(dirs["bmp_edge"], "odd_width_32bit.bmp"), odd, bits=32)
    write_ppm(os.path.join(dirs["ppm"], "odd_width_four_colors.ppm"), odd)

    print(f"Wrote test images to {root}/")


if __name__ == "__main__":
    main()
