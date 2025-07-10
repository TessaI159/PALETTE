#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "Color.h"

#define COLOR_SPACE_BIT(space) (1u << (space))
#define square(x) ((x) * (x))
#define second_sursolid(x) ((x) * (x) * (x) * (x) * (x) * (x) * (x))
#define DEG2RAD(x) (((x) / (180.0f)) * (M_PI))
#define RAD2DEG(x) (((x) / M_PI) * (180.0f))

/* TODO (Tess): Remove any unnecessary indirections */

/* Constants for a standard D65/2 illuminant */
static const double X2	  = 0.95047;
static const double Y2	  = 1.000;
static const double Z2	  = 1.08883;
static const double DELTA = (6.0 / 29.0) * (6.0 / 29.0) * (6.0 / 29.0);
#ifndef M_PI
static const double M_PI = 3.14159265358979323846;
#endif

/* Mark/check the validity of colors */
static inline void Color_mark_space(struct Color   *color,
				    enum ColorSpace space) {
	color->valid_spaces |= COLOR_SPACE_BIT(space);
}

#ifdef PALETTE_DEBUG
int Color_has_space(const struct Color *color, enum ColorSpace space) {
#else
static inline int Color_has_space(const struct Color *color,
				  enum ColorSpace     space) {
#endif
	return (color->valid_spaces & COLOR_SPACE_BIT(space)) != 0;
}

static inline float linearize(float channel) {
	return (channel > 0.04045f) ? powf((channel + 0.055f) / 1.055f, 2.4f)
				    : channel / 12.92f;
}

static inline float delinearize(float channel) {
	return (channel > 0.0031308f)
		   ? 1.055f * powf(channel, 1.0f / 2.4f) - 0.055f
		   : 12.92f * channel;
}

struct Color Color_create(const float r, const float g, const float b) {
	struct Color color = {0};
	color.srgb.r	   = linearize(r / 255.0f);
	color.srgb.g	   = linearize(g / 255.0f);
	color.srgb.b	   = linearize(b / 255.0f);
	Color_mark_space(&color, COLOR_SRGB);
	return color;
}

struct Color Color_create_norm(const float r, const float g, const float b) {
	struct Color color;
	color.srgb.r	   = linearize(r);
	color.srgb.g	   = linearize(g);
	color.srgb.b	   = linearize(b);
	color.valid_spaces = 0;
	Color_mark_space(&color, COLOR_SRGB);
	return color;
}

static void convert_srgb_to_cielab(struct Color *color) {
	const struct sRGB srgb = color->srgb;

	double r = srgb.r;
	double g = srgb.g;
	double b = srgb.b;

	double x = fma(r, 0.4124564, fma(g, 0.3575761, fma(b, 0.1804375, 0.0)));
	double y = fma(r, 0.2126729, fma(g, 0.7151522, fma(b, 0.0721750, 0.0)));
	double z = fma(r, 0.0193339, fma(g, 0.1191920, fma(b, 0.9503041, 0.0)));

	x /= X2;
	y /= Y2;
	z /= Z2;

	double fx = x > DELTA ? cbrt(x) : fma(7.787, x, 16.0 / 116.0);
	double fy = y > DELTA ? cbrt(y) : fma(7.787, y, 16.0 / 116.0);
	double fz = z > DELTA ? cbrt(z) : fma(7.787, z, 16.0 / 116.0);

	color->cielab.l = fma(116.0, fy, -16.0);
	color->cielab.a = 500.0 * (fx - fy);
	color->cielab.b = 200.0 * (fy - fz);

	Color_mark_space(color, COLOR_CIELAB);
}

static void convert_srgb_to_oklab(struct Color *color) {
	const struct sRGB srgb = color->srgb;

	double l =
	    fma(0.4122214708, srgb.r,
		fma(0.5363325363, srgb.g, fma(0.0514459929, srgb.b, 0.0)));
	double m =
	    fma(0.2119034982, srgb.r,
		fma(0.6806995451, srgb.g, fma(0.1073969566, srgb.b, 0.0)));
	double s =
	    fma(0.0883024619, srgb.r,
		fma(0.2817188376, srgb.g, fma(0.6299787005, srgb.b, 0.0)));

	l = cbrtf(l);
	m = cbrtf(m);
	s = cbrtf(s);

	color->oklab.l = fma(0.2104542553, l,
			     fma(0.7936177850, m, fma(-0.0040720468, s, 0.0)));
	color->oklab.a = fma(1.9779984951, l,
			     fma(-2.4285922050, m, fma(0.4505937099, s, 0.0)));
	color->oklab.b = fma(0.0259040371, l,
			     fma(0.7827717662, m, fma(-0.8086757660, s, 0.0)));

	Color_mark_space(color, COLOR_OKLAB);
}
#ifdef PALETTE_DEBUG
void convert_srgb_to_grayscale(struct Color *color) {
#else
static inline void convert_srgb_to_grayscale(struct Color *color) {
#endif
	struct sRGB srgb = color->srgb;

	float l =
	    fma(0.2126, srgb.r, fma(0.7152, srgb.g, fma(0.0722, srgb.b, 0)));

	color->grayscale.l = l;
	Color_mark_space(color, COLOR_GRAY);
}

/* Color difference functions begin here */
static inline float delta_ok_diff_fast(const struct okLAB *sam,
				       const struct okLAB *ref) {
	return sqrtf((sam->l - ref->l) * (sam->l - ref->l) +
		     (sam->a - ref->a) * (sam->a - ref->a) +
		     (sam->b - ref->b) * (sam->b - ref->b));
}

static inline float delta_cie76_diff_fast(const struct cieLAB *sam,
					  const struct cieLAB *ref) {
	return sqrtf(fmaf(sam->l - ref->l, sam->l - ref->l,
			  fmaf(sam->a - ref->a, sam->a - ref->a,
			       fmaf(sam->b - ref->b, sam->b - ref->b, 0))));
}

static inline float delta_cie94_diff_fast(const struct cieLAB *sam,
					  const struct cieLAB *ref) {
	const double K1 = 0.045;
	const double K2 = 0.015;

	const double dl = sam->l - ref->l;
	const double da = sam->a - ref->a;
	const double db = sam->b - ref->b;

	const double C1 = sqrt(square(sam->a) + square(sam->b));
	const double C2 = sqrt(square(ref->a) + square(ref->b));
	const double dC = C1 - C2;

	const double dH2 = square(da) + square(db) - square(dC);

	const double sl = 1.0;
	const double sc = fmaf(K1, C1, 1.0);
	const double sh = fmaf(K2, C1, 1.0);

	const double term_L = dl / sl;
	const double term_C = dC / sc;
	const double term_H = sqrt(dH2) / sh;

	return sqrt(
	    fma(term_L, term_L, fma(term_C, term_C, fma(term_H, term_H, 0.0))));
}

static inline float delta_ciede2000_diff_fast(const struct cieLAB *sam,
					      const struct cieLAB *ref) {

	const double kl = 1;
	const double kh = 1;
	const double kc = 1;
	double	     l1 = sam->l, a1 = sam->a, b1 = sam->b;
	double	     l2 = ref->l, a2 = ref->a, b2 = ref->b;

	/* Step 1 */
	double c1 = sqrtf(square(a1) + square(b1));
	double c2 = sqrtf(square(a2) + square(b2));

	double avgc = (c1 + c2) / 2.0f;
	/* Saves calculating this twice */
	double g =
	    0.5f * (1.0f - sqrtf(second_sursolid(avgc) /
				 (second_sursolid(avgc) + 6103515625.0)));
	double a1p = (1.0 + g) * a1;
	double a2p = (1.0 + g) * a2;
	double c1p = sqrt(square(a1p) + square(b1));
	double c2p = sqrt(square(a2p) + square(b2));
	double h1p = (b1 == 0.0 && a1p == 0.0)
			 ? 0
			 : fmod(RAD2DEG(atan2(b1, a1p)) + 360.0, 360.0);
	double h2p = (b2 == 0.0f && a2p == 0.0f)
			 ? 0
			 : fmod(RAD2DEG(atan2(b2, a2p)) + 360.0, 360.0);

	/* Step 2 */
	double dlp = l2 - l1;
	double dcp = c2p - c1p;

	double dhp;
	if (c1p == 0.0 || c2p == 0.0) {
		dhp = 0.0f;
	} else {
		if (fabs(h2p - h1p) <= 180.0) {
			dhp = h2p - h1p;
		} else if (h2p - h1p > 180.0) {
			dhp = h2p - h1p - 360.0;
		} else if (h2p - h1p < -180.0) {
			dhp = h2p - h1p + 360.0;
		} else {
			return 0.0;
		}
	}

	double dHp = 2.0f * sqrt(c1p * c2p) * sin(DEG2RAD(dhp / 2.0));

	/* Step 3 */
	double avglp = (l1 + l2) / 2.0;
	double avgcp = (c1p + c2p) / 2.0;

	double avghp;

	if (c1p == 0.0 || c2p == 0.0) {
		avghp = h1p + h2p;
	} else {
		if (fabs(h1p - h2p) <= 180.0) {
			avghp = (h1p + h2p) / 2.0;
		} else {
			if (h1p + h2p < 360.0) {
				avghp = (h1p + h2p + 360.0) / 2.0;
			} else {
				avghp = (h1p + h2p - 360.0) / 2.0;
			}
		}
	}

	double t = 1.0 - 0.17 * cos(DEG2RAD(avghp - 30.0)) +
		   0.24 * cos(DEG2RAD(2.0 * avghp)) +
		   0.32 * cos(DEG2RAD(3.0 * avghp + 6.0)) -
		   0.20 * cos(DEG2RAD(4.0 * avghp - 63.0));
	double dtheta = 30.0 * exp(-(square((avghp - 275.0) / 25.0)));
	double rc     = 2.0 * sqrt(second_sursolid(avgcp) /
				   (second_sursolid(avgcp) + 6103515625.0f));
	double sl     = 1.0 + ((0.015 * (square(avglp - 50.0))) /
			       sqrt(20.0 + (square(avglp - 50.0))));
	double sc     = 1.0 + 0.045 * avgcp;
	double sh     = 1.0 + 0.015 * avgcp * t;
	double rt     = -sin(DEG2RAD(2.0 * dtheta)) * rc;
	double terml  = dlp / (kl * sl);
	double termc  = dcp / (kc * sc);
	double termh  = dHp / (kh * sh);

	return sqrt(square(terml) + square(termc) + square(termh) +
		    rt * termc * termh);
}

float delta_ok_diff(struct Color *sam, struct Color *ref) {
	if (!Color_has_space(sam, COLOR_OKLAB)) {
		convert_srgb_to_oklab(sam);
	}
	if (!Color_has_space(ref, COLOR_OKLAB)) {
		convert_srgb_to_oklab(ref);
	}
	return delta_ok_diff_fast(&sam->oklab, &ref->oklab);
}

float delta_cie76_diff(struct Color *sam, struct Color *ref) {
	if (!Color_has_space(sam, COLOR_CIELAB)) {
		convert_srgb_to_cielab(sam);
	}
	if (!Color_has_space(ref, COLOR_CIELAB)) {
		convert_srgb_to_cielab(ref);
	}
	return delta_cie76_diff_fast(&sam->cielab, &ref->cielab);
}

float delta_cie94_diff(struct Color *sam, struct Color *ref) {
	if (!Color_has_space(sam, COLOR_CIELAB)) {
		convert_srgb_to_cielab(sam);
	}
	if (!Color_has_space(ref, COLOR_CIELAB)) {
		convert_srgb_to_cielab(ref);
	}

	return delta_cie94_diff_fast(&sam->cielab, &ref->cielab);
}

float delta_ciede2000_diff(struct Color *sam, struct Color *ref) {
	if (!Color_has_space(sam, COLOR_CIELAB)) {
		convert_srgb_to_cielab(sam);
	}
	if (!Color_has_space(ref, COLOR_CIELAB)) {
		convert_srgb_to_cielab(ref);
	}
	return delta_ciede2000_diff_fast(&sam->cielab, &ref->cielab);
}

void Color_calc_spaces(struct Color *color) {
	convert_srgb_to_oklab(color);
	convert_srgb_to_cielab(color);
	convert_srgb_to_grayscale(color);
}

#ifdef PALETTE_DEBUG
void Color_print(struct Color *color) {
	printf("Linear sRGB: (%f, %f, %f)\n", color->srgb.r, color->srgb.g,
	       color->srgb.b);
	if (!Color_has_space(color, COLOR_OKLAB)) {
		convert_srgb_to_oklab(color);
	}
	printf("okLAB: (%f, %f, %f)\n", color->oklab.l, color->oklab.a,
	       color->oklab.b);
	if (!Color_has_space(color, COLOR_CIELAB)) {
		convert_srgb_to_cielab(color);
	}
	printf("cieLAB: (%f, %f, %f)\n", color->cielab.l, color->cielab.a,
	       color->cielab.b);
	if (!Color_has_space(color, COLOR_GRAY)) {
		convert_srgb_to_grayscale(color);
	}
	printf("Grayscale: %f\n", color->grayscale.l);
}
#endif
