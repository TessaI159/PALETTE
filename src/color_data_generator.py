import csv
import argparse
import numpy as np
from colour import XYZ_to_Oklab
from colour.difference import (
    delta_E_CIE1976,
    delta_E_CIE1994,
    delta_E_CIE2000,
)

# Reference ColorChecker sRGB values under D65
reference_colors_srgb = [
    (0.451, 0.322, 0.267), (0.761, 0.588, 0.510), (0.384, 0.478, 0.616),
    (0.341, 0.424, 0.263), (0.522, 0.502, 0.694), (0.404, 0.741, 0.667),
    (0.839, 0.494, 0.173), (0.314, 0.357, 0.651), (0.757, 0.353, 0.388),
    (0.369, 0.235, 0.424), (0.616, 0.737, 0.251), (0.878, 0.639, 0.180),
    (0.220, 0.239, 0.588), (0.275, 0.580, 0.286), (0.686, 0.212, 0.235),
    (0.906, 0.780, 0.122), (0.733, 0.337, 0.584), (0.031, 0.522, 0.631),
    (0.953, 0.953, 0.949), (0.784, 0.784, 0.784), (0.627, 0.627, 0.627),
    (0.478, 0.478, 0.475)
]

# D65 white point
Xn, Yn, Zn = 0.95047, 1.00000, 1.08883
delta = (6.0 / 29.0) ** 3

# RGB → XYZ matrix for linear RGB
RGB_TO_XYZ_MATRIX = np.array([
    [0.4124564, 0.3575761, 0.1804375],
    [0.2126729, 0.7151522, 0.0721750],
    [0.0193339, 0.1191920, 0.9503041]
])

def linearize_srgb(srgb):
    return np.where(srgb <= 0.04045,
                    srgb / 12.92,
                    ((srgb + 0.055) / 1.055) ** 2.4)

def f(t):
    return np.where(t > delta,
                    np.cbrt(t),
                    7.787 * t + 16.0 / 116.0)

def linear_rgb_to_xyz(rgb_lin):
    return RGB_TO_XYZ_MATRIX @ rgb_lin

def xyz_to_lab(xyz):
    x, y, z = xyz[0] / Xn, xyz[1] / Yn, xyz[2] / Zn
    fx, fy, fz = f(np.array([x, y, z]))
    L = 116.0 * fy - 16.0
    a = 500.0 * (fx - fy)
    b = 200.0 * (fy - fz)
    return np.array([L, a, b])

def euclidean_diff(c1, c2):
    return np.linalg.norm(c1 - c2)

def generate_color_data():
    srgb = np.array(reference_colors_srgb)
    linear_srgb = linearize_srgb(srgb)

    xyz = np.array([linear_rgb_to_xyz(c) for c in linear_srgb])
    lab = np.array([xyz_to_lab(x) for x in xyz])
    oklab = np.array([XYZ_to_Oklab(x) for x in xyz])

    grayscale = np.dot(linear_srgb, [0.2126, 0.7152, 0.0722])

    return srgb, linear_srgb, lab, oklab, grayscale

def write_reference_csv(path, srgb, linear_srgb, lab, oklab, grayscale):
    fieldnames = [
        "index",
        "srgb_r", "srgb_g", "srgb_b",
        "linear_r", "linear_g", "linear_b",
        "lab_l", "lab_a", "lab_b",
        "oklab_l", "oklab_a", "oklab_b",
        "grayscale"
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i in range(len(srgb)):
            writer.writerow({
                "index": i,
                "srgb_r": srgb[i][0], "srgb_g": srgb[i][1], "srgb_b": srgb[i][2],
                "linear_r": linear_srgb[i][0], "linear_g": linear_srgb[i][1], "linear_b": linear_srgb[i][2],
                "lab_l": lab[i][0], "lab_a": lab[i][1], "lab_b": lab[i][2],
                "oklab_l": oklab[i][0], "oklab_a": oklab[i][1], "oklab_b": oklab[i][2],
                "grayscale": grayscale[i]
            })

def write_diff_csv(path, lab, oklab):
    fieldnames = [
        "color1_index", "color2_index",
        "oklab_diff", "cie76_diff", "cie94_diff", "ciede2000_diff"
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i in range(len(lab)):
            for j in range(i + 1, len(lab)):
                writer.writerow({
                    "color1_index": i,
                    "color2_index": j,
                    "oklab_diff": euclidean_diff(oklab[i], oklab[j]),
                    "cie76_diff": delta_E_CIE1976(lab[i], lab[j]),
                    "cie94_diff": delta_E_CIE1994(lab[i], lab[j]),
                    "ciede2000_diff": delta_E_CIE2000(lab[i], lab[j]),
                })

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", type=str, default="color_output")
    args = parser.parse_args()

    srgb, linear_srgb, lab, oklab, grayscale = generate_color_data()

    ref_path = f"{args.output_prefix}_reference_colors.csv"
    diff_path = f"{args.output_prefix}_pairwise_differences.csv"

    write_reference_csv(ref_path, srgb, linear_srgb, lab, oklab, grayscale)
    write_diff_csv(diff_path, lab, oklab)

    print(f"[✓] Reference colors written to: {ref_path}")
    print(f"[✓] Pairwise differences written to: {diff_path}")

if __name__ == "__main__":
    main()
