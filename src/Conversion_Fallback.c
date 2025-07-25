#include <math.h>

#include "Color.h"

#define CUBE(x) ((x) * (x) * (x))

/* Standard D65 2 degree illuminant constants */
#define X2 95.047
#define Y2 100.0
#define Z2 108.883

/* Using doubles internally and then casting to float when assigning back to
   struct Color adds just a touch of precision due to rounding, at the cost of
   <= 24 bytes */

static inline void srgb_to_xyz(double *restrict r, double *restrict g,
			       double *restrict b, double *restrict x,
			       double *restrict y, double *restrict z) {
	*r *= 100.0f;
	*g *= 100.0f;
	*b *= 100.0f;
	*x = *r * 0.4124 + *g * 0.3576 + *b * 0.1805;
	*y = *r * 0.2126 + *g * 0.7152 + *b * 0.0722;
	*z = *r * 0.0193 + *g * 0.1192 + *b * 0.9505;
}

static inline void xyz_to_srgb(double *restrict x, double *restrict y,
			       double *restrict z, double *restrict r,
			       double *restrict g, double *restrict b) {
	*r = *x * 3.2406 + *y * -1.5372 + *z * -0.4986;
	*g = *x * -0.9689 + *y * 1.8758 + *z * 0.0415;
	*b = *x * 0.0557 + *y * -0.2040 + *z * 1.0570;
}

/* Assumes linear srgb input */
void srgb_to_oklab_fb(struct Color *restrict colors) {
	const uint32_t num_col = colors->num;
	for (uint32_t i = 0; i < num_col; ++i) {
		double r = colors->alpha[i];
		double g = colors->beta[i];
		double b = colors->gamma[i];
		double l =
		    0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b;
		double m =
		    0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b;
		double s =
		    0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b;
		l = cbrt(l);
		m = cbrt(m);
		s = cbrt(s);

		colors->alpha[i] =
		    0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s;
		colors->beta[i] =
		    1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s;
		colors->gamma[i] =
		    0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s;
	}
}

/* Assumes linear srgb input */
void srgb_to_cielab_fb(struct Color *restrict colors) {
	const uint32_t num_col = colors->num;
	for (uint32_t i = 0; i < num_col; ++i) {
		const double r = colors->alpha[i] * 100.0;
		const double g = colors->beta[i] * 100.0;
		const double b = colors->gamma[i] * 100.0;
		double	     x, y, z;
		srgb_to_xyz(r, g, b, &x, &y, &z);
		x /= X2;
		y /= Y2;
		z /= Z2;

		x = x > 0.008856 ? cbrt(x) : (7.787 * x) + (16.0 / 116.0);
		y = y > 0.008856 ? cbrt(y) : (7.787 * y) + (16.0 / 116.0);
		z = z > 0.008856 ? cbrt(z) : (7.787 * z) + (16.0 / 116.0);

		colors->alpha[i] = (116.0 * y) - 16.0;
		colors->beta[i]	 = 500.0 * (x - y);
		colors->gamma[i] = 200.0 * (y - z);
	}
}

/* Provides linear srgb output */

void oklab_to_srgb_fb(struct Color *restrict colors) {
	const uint32_t num_col = colors->num;
	for (uint32_t i = 0; i < num_col; ++i) {
		const double l	= colors->alpha[i];
		const double a	= colors->beta[i];
		const double b	= colors->gamma[i];
		double	     l_ = l + 0.3963377774 * a + 0.2158037573 * b;
		double	     m_ = l - 0.1055613458 * a - 0.0638541728 * b;
		double	     s_ = l - 0.0894841775 * a - 1.2914855480 * b;

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

/* Provides linear srgb output */
void cielab_to_srgb_fb(struct Color *restrict colors) {
	const uint32_t num_col = colors->num;
	for (uint32_t i = 0; i < num_col; ++i) {
		const double l = colors->alpha[i];
		const double a = colors->beta[i];
		const double b = colors->gamma[i];
		double	     y = (l + 16.0) / 116.0;
		double	     x = a / 500.0 + y;
		double	     z = y - b / 200.0;

		x = CUBE(x) > 0.008856 ? CUBE(x) : (x - 16.0 / 116.0) / 7.787;
		y = CUBE(y) > 0.008856 ? CUBE(y) : (y - 16.0 / 116.0) / 7.787;
		z = CUBE(z) > 0.008856 ? CUBE(z) : (z - 16.0 / 116.0) / 7.787;

		colors->alpha[i] = x * X2;
		colors->beta[i]	 = y * Y2;
		colors->gamma[i] = z * Z2;
	}
}
