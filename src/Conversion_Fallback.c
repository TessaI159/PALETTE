#include <math.h>
#include <stdint.h>

#include "Color.h"
#include "Conversion.h"

#define CUBE(x) ((x) * (x) * (x))
#define SQUARE(x) ((x) * (x))

/* Standard D65 2 degree illuminant constants */
#define INV_X2 (1.0 / 95.047)
#define INV_Y2 (1.0 / 100.0)
#define INV_Z2 (1.0 / 108.883)

#define X2 95.047
#define Y2 100.0
#define Z2 108.883

void srgb_linearize(struct Color *restrict colors) {
  const static double threshold = 0.04045;
  const uint32_t num_col = colors->num;
  for (uint32_t i = 0; i < num_col; ++i) {

    const double r = colors->alpha[i];
    const double g = colors->beta[i];
    const double b = colors->gamma[i];
    colors->alpha[i] =
        (r > threshold) ? pow((r + 0.055) / 1.055, 2.4) : r / 12.92;
    colors->beta[i] =
        (g > threshold) ? pow((g + 0.055) / 1.055, 2.4) : g / 12.92;
    colors->gamma[i] =
        (b > threshold) ? pow((b + 0.055) / 1.055, 2.4) : b / 12.92;
  }
}

void srgb_gamma_correct(struct Color *restrict const colors) {
  const static double threshold = 0.0031308;
  const uint32_t num_col = colors->num;
  for (uint32_t i = 0; i < num_col; ++i) {
    const double r = colors->alpha[i];
    const double g = colors->beta[i];
    const double b = colors->gamma[i];

    colors->alpha[i] =
        (r > threshold) ? 1.055 * pow(r, 1.0 / 2.4) - 0.055 : 12.92 * r;
    colors->beta[i] =
        (g > threshold) ? 1.055 * pow(g, 1.0 / 2.4) - 0.055 : 12.92 * g;
    colors->gamma[i] =
        (b > threshold) ? 1.055 * pow(b, 1.0 / 2.4) - 0.055 : 12.92 * b;
  }
}

static inline void srgb_to_xyz(struct Color *restrict const colors) {
  const uint32_t num_col = colors->num;
  for (uint32_t i = 0; i < num_col; ++i) {
    const double r = colors->alpha[i] * 100.0;
    const double g = colors->beta[i] * 100.0;
    const double b = colors->gamma[i] * 100.0;
    const double x = r * 0.4124564 + g * 0.3575761 + b * 0.1804375;
    const double y = r * 0.2126729 + g * 0.7151522 + b * 0.0721750;
    const double z = r * 0.0193339 + g * 0.1191920 + b * 0.9503041;
    colors->alpha[i] = x;
    colors->beta[i] = y;
    colors->gamma[i] = z;
  }
}

static inline void xyz_to_srgb(struct Color *restrict const colors) {
  const uint32_t num_col = colors->num;
  for (uint32_t i = 0; i < num_col; ++i) {
    const double x = colors->alpha[i] * 0.01;
    const double y = colors->beta[i] * 0.01;
    const double z = colors->gamma[i] * 0.01;
    const double r = x * 3.2404542 + y * -1.5371385 + z * -0.4985314;
    const double g = x * -0.9692660 + y * 1.8760108 + z * 0.0415560;
    const double b = x * 0.0556434 + y * -0.2040259 + z * 1.0572252;
    colors->alpha[i] = r;
    colors->beta[i] = g;
    colors->gamma[i] = b;
  }
}

/* Linear sRGB [0,1] -> Oklab */
void srgb_to_oklab_fb(struct Color *restrict colors) {
  const uint32_t num_col = colors->num;
  for (uint32_t i = 0; i < num_col; ++i) {
    const double r = colors->alpha[i];
    const double g = colors->beta[i];
    const double b = colors->gamma[i];
    double l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b;
    double m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b;
    double s = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b;
    l = cbrt(l);
    m = cbrt(m);
    s = cbrt(s);

    colors->alpha[i] = 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s;
    colors->beta[i] = 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s;
    colors->gamma[i] = 0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s;
  }
}

/* Oklab -> Linear sRGB [0, 1]*/
void oklab_to_srgb_fb(struct Color *restrict colors) {
  const uint32_t num_col = colors->num;
  for (uint32_t i = 0; i < num_col; ++i) {
    const double l = colors->alpha[i];
    const double a = colors->beta[i];
    const double b = colors->gamma[i];
    double l_ = l + 0.3963377774 * a + 0.2158037573 * b;
    double m_ = l - 0.1055613458 * a - 0.0638541728 * b;
    double s_ = l - 0.0894841775 * a - 1.2914855480 * b;

    l_ = CUBE(l_);
    m_ = CUBE(m_);
    s_ = CUBE(s_);

    colors->alpha[i] =
        4.0767416621 * l_ - 3.3077115913 * m_ + 0.2309699292 * s_;
    colors->beta[i] =
        -1.2684380046 * l_ + 2.6097574011 * m_ - 0.3413193965 * s_;
    colors->gamma[i] =
        -0.0041960863 * l_ - 0.7034186147 * m_ + 1.7076147010 * s_;
  }
}

/* Linear sRGB [0, 1] -> CIELAB */
void srgb_to_cielab_fb(struct Color *restrict colors) {
  const uint32_t num_col = colors->num;
  srgb_to_xyz(colors);
  const static double delta = 6.0 / 29.0;
  const static double delta2 = SQUARE(delta);
  const static double delta3 = CUBE(delta);

  for (uint32_t i = 0; i < num_col; ++i) {
    double x = colors->alpha[i] * INV_X2;
    double y = colors->beta[i] * INV_Y2;
    double z = colors->gamma[i] * INV_Z2;

    x = (x > delta3) ? cbrt(x) : (x / (3.0 * delta2) + 4.0 / 29.0);
    y = (y > delta3) ? cbrt(y) : (y / (3.0 * delta2) + 4.0 / 29.0);
    z = (z > delta3) ? cbrt(z) : (z / (3.0 * delta2) + 4.0 / 29.0);

    colors->alpha[i] = (116.0 * y) - 16.0;
    colors->beta[i] = 500.0 * (x - y);
    colors->gamma[i] = 200.0 * (y - z);
  }
}

/* CIELAB -> Linear sRGB [0, 1] */
void cielab_to_srgb_fb(struct Color *restrict colors) {
  const uint32_t num_col = colors->num;
  const static double delta = 6.0 / 29.0;
  const static double delta2 = SQUARE(delta);
  const static double delta3 = CUBE(delta);
  for (uint32_t i = 0; i < num_col; ++i) {
    const double l = colors->alpha[i];
    const double a = colors->beta[i];
    const double b = colors->gamma[i];

    const double fy = (l + 16.0) / 116.0;
    const double fx = a / 500.0 + fy;
    const double fz = fy - b / 200.0;

    const double x =
        (CUBE(fx) > delta3) ? CUBE(fx) : (fx - 16.0 / 116.0) * 3.0 * delta2;
    const double y =
        (CUBE(fy) > delta3) ? CUBE(fy) : (fy - 16.0 / 116.0) * 3.0 * delta2;
    const double z =
        (CUBE(fz) > delta3) ? CUBE(fz) : (fz - 16.0 / 116.0) * 3.0 * delta2;

    colors->alpha[i] = x * X2;
    colors->beta[i] = y * Y2;
    colors->gamma[i] = z * Z2;
  }
  xyz_to_srgb(colors);
}
