#include <math.h>

#include "Color.h"

#define CUBE(x) ((x) * (x) * (x))
/* Micro-optimizing out division */
#define INV_255 0.00392156862745098
#define INV_1_055 0.9478672985781991
#define INV_12_92 0.07739938080495357

/* Standard D65 2 degree illuminant constants */
#define INV_X2 1.0521110608435826
#define INV_Y2 1.0
#define INV_Z2 .9184170164304804

#define X2 95.047
#define Y2 100.0
#define Z2 108.883

/* Using doubles internally and then casting to float when assigning back to
   struct Color adds just a touch of precision due to rounding, at the cost of
   <= 72 bytes */
static inline void linearize(struct Color *restrict const colors) {
	const double threshold = 0.04045;
	uint16_t     num_col   = colors->num;
	for (uint16_t i = 0; i < num_col; ++i) {
		colors->alpha[i] *= INV_255;
		colors->beta[i] *= INV_255;
		colors->gamma[i] *= INV_255;

		double tmpr = colors->alpha[i];
		double tmpg = colors->beta[i];
		double tmpb = colors->gamma[i];

		colors->alpha[i] = tmpr > threshold
				       ? pow(((tmpr + 0.055) * INV_1_055), 2.4)
				       : tmpr * INV_12_92;
		colors->beta[i]	 = tmpg > threshold
				       ? pow(((tmpg + 0.055) * INV_1_055), 2.4)
				       : tmpg * INV_12_92;
		colors->gamma[i] = tmpb > threshold
				       ? pow(((tmpb + 0.055) * INV_1_055), 2.4)
				       : tmpb * INV_12_92;
	}
}

static inline void gamma_correct(struct Color *restrict const colors) {
	const double threshold = 0.0031308;
	uint16_t     num_col   = colors->num;
	for (uint16_t i = 0; i < num_col; ++i) {
		double tmpr = colors->alpha[i];
		double tmpg = colors->beta[i];
		double tmpb = colors->gamma[i];

		colors->alpha[i] =
		    tmpr > threshold
			? 1.055 * pow(tmpr, 0.4166666666666663) - 0.055
			: 12.92 * tmpr;
		colors->beta[i] =
		    tmpg > threshold
			? 1.055 * pow(tmpg, 0.4166666666666663) - 0.055
			: 12.92 * tmpg;
		colors->gamma[i] =
		    tmpb > threshold
			? 1.055 * pow(tmpb, 0.4166666666666663) - 0.055
			: 12.92 * tmpb;
	}
}

static inline void srgb_to_xyz(struct Color *restrict const colors) {
	uint16_t num_col = colors->num;
	linearize(colors);
	for (uint16_t i = 0; i < num_col; ++i) {
		colors->alpha[i] *= 100.0;
		colors->beta[i] *= 100.0;
		colors->gamma[i] *= 100.0;
		float tmpx = colors->alpha[i] * 0.4124 +
			     colors->beta[i] * 0.3576 +
			     colors->gamma[i] * 0.1805;
		float tmpy = colors->alpha[i] * 0.2126 +
			     colors->beta[i] * 0.7152 +
			     colors->gamma[i] * 0.0722;
		float tmpz = colors->alpha[i] * 0.0193 +
			     colors->beta[i] * 0.1192 +
			     colors->gamma[i] * 0.9505;
		colors->alpha[i] = tmpx;
		colors->beta[i]	 = tmpy;
		colors->gamma[i] = tmpz;
	}
}

static inline void xyz_to_srgb(struct Color *restrict const colors) {
	uint16_t num_col = colors->num;
	for (uint16_t i = 0; i < num_col; ++i) {
		colors->alpha[i] *= 0.01;
		colors->beta[i] *= 0.01;
		colors->gamma[i] *= 0.01;
		float tmpr = colors->alpha[i] * 3.2406 +
			     colors->beta[i] * -1.5372 +
			     colors->gamma[i] * -0.4986;
		float tmpg = colors->alpha[i] * -0.9689 +
			     colors->beta[i] * 1.8758 +
			     colors->gamma[i] * 0.0415;
		float tmpb = colors->alpha[i] * 0.0557 +
			     colors->beta[i] * -0.2040 +
			     colors->gamma[i] * 1.0570;
		colors->alpha[i] = tmpr;
		colors->beta[i]	 = tmpg;
		colors->gamma[i] = tmpb;
	}
}

void srgb_to_oklab_fb(struct Color *restrict colors) {
	const uint32_t num_col = colors->num;
	for (uint32_t i = 0; i < num_col; ++i) {
		const double r = colors->alpha[i];
		const double g = colors->beta[i];
		const double b = colors->gamma[i];
		double	     l =
		    0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b;
		double m =
		    0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b;
		double s =
		    0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b;
		l = cbrt(l);
		m = cbrt(m);
		s = cbrt(s);
		/* cbrt has domain [0, 1] */

		colors->alpha[i] =
		    0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s;
		colors->beta[i] =
		    1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s;
		colors->gamma[i] =
		    0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s;
	}
}

void srgb_to_cielab_fb(struct Color *restrict colors) {
	const uint32_t num_col = colors->num;
	srgb_to_xyz(colors);
	for (uint32_t i = 0; i < num_col; ++i) {
		double x = colors->alpha[i];
		double y = colors->beta[i];
		double z = colors->gamma[i];

		x *= INV_X2;
		y *= INV_Y2;
		z *= INV_Z2;

		x = x > 0.008856 ? cbrt(x)
				 : (7.787 * x) + (0.13793103448275862);
		y = y > 0.008856 ? cbrt(y)
				 : (7.787 * y) + (0.13793103448275862);
		z = z > 0.008856 ? cbrt(z)
				 : (7.787 * z) + (0.13793103448275862);

		colors->alpha[i] = (116.0 * y) - 16.0;
		colors->beta[i]	 = 500.0 * (x - y);
		colors->gamma[i] = 200.0 * (y - z);
	}
}

/* Provides linear srgb output */
void oklab_to_srgb_fb(struct Color *restrict colors) {
	uint32_t num_col = colors->num;
	for (uint32_t i = 0; i < num_col; ++i) {
		double l  = colors->alpha[i];
		double a  = colors->beta[i];
		double b  = colors->gamma[i];
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

/* Provides linear srgb output */
void cielab_to_srgb_fb(struct Color *restrict colors) {
	uint32_t num_col = colors->num;
	for (uint32_t i = 0; i < num_col; ++i) {
		const double l	= colors->alpha[i];
		const double a	= colors->beta[i];
		const double b	= colors->gamma[i];
		double	     y	= (l + 16.0) / 116.0;
		double	     x	= a / 500.0 + y;
		double	     z	= y - b / 200.0;
		const double x3 = CUBE(x);
		const double y3 = CUBE(y);
		const double z3 = CUBE(z);

		x = x3 > 0.008856 ? x3 : (x - 0.13793103448275862) / 7.787;
		y = y3 > 0.008856 ? y3 : (y - 0.13793103448275862) / 7.787;
		z = z3 > 0.008856 ? z3 : (z - 0.13793103448275862) / 7.787;

		colors->alpha[i] = x * X2;
		colors->beta[i]	 = y * Y2;
		colors->gamma[i] = z * Z2;
	}
	xyz_to_srgb(colors);
}
