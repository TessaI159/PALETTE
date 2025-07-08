#define _USE_MATH_DEFINES
#include <math.h>

#include "Color.h"

#define COLOR_SPACE_BIT(space) (1u << (space))
#define square(x) ((x)*(x))

/* Constants for a standard D65/2 illuminant */
static const double X2	  = 0.95047;
static const double Y2	  = 1.000;
static const double Z2	  = 1.08883;
static const double DELTA = (6.0 / 29.0) * (6.0 / 29.0) * (6.0 / 29.0);
#ifndef M_PI
static const double M_PI = 3.14159265358979323846;
#endif

/* Mark/check the validity of colors */
/* NOTE (Tess): Currently unused due to lack of supported color spaces*/
static inline void Color_mark_space(struct Color   *color,
				    enum ColorSpace space) {
	color->valid_spaces |= COLOR_SPACE_BIT(space);
}

static inline int Color_has_space(const struct Color color,
				  enum ColorSpace    space) {
	return (color.valid_spaces & COLOR_SPACE_BIT(space)) != 0;
}

static inline double linearize(double channel) {
	return (channel > 0.04045) ? pow((channel + 0.055) / 1.055, 2.4)
				   : channel / 12.92;
}

static inline double delinearize(double channel) {
	return (channel > 0.0031308) ? 1.055 * pow(channel, 1.0 / 2.4) - 0.055
				     : 12.92 * channel;
}

struct Color Color_create(const double r, const double g, const double b) {
	struct Color color;
	color.srgb.r = linearize(r / 255.0);
	color.srgb.g = linearize(g / 255.0);
	color.srgb.b = linearize(b / 255.0);
	Color_mark_space(&color, COLOR_SRGB);
	return color;
}

struct Color Color_create_norm(const double r, const double g, const double b) {
	struct Color color;
	color.srgb.r = linearize(r);
	color.srgb.g = linearize(g);
	color.srgb.b = linearize(b);
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

	l = cbrt(l);
	m = cbrt(m);
	s = cbrt(s);

	color->oklab.l = fma(0.2104542553, l,
			     fma(0.7936177850, m, fma(-0.0040720468f, s, 0.0)));
	color->oklab.a = fma(1.9779984951, l,
			     fma(-2.4285922050, m, fma(0.4505937099, s, 0.0)));
	color->oklab.b = fma(0.0259040371, l,
			     fma(0.7827717662, m, fma(-0.8086757660, s, 0.0)));

	Color_mark_space(color, COLOR_OKLAB);
}

static inline void convert_srgb_to_grayscale(struct Color *color) {
	struct sRGB srgb = color->srgb;

	double l =
	    fma(0.2126, srgb.r, fma(0.7152, srgb.g, fma(0.0722, srgb.b, 0)));

	color->grayscale.l = l;
}

/* Color difference functions begin here */
static inline double delta_ok_diff_fast(const struct okLAB sam,
					const struct okLAB ref) {
	return sqrt((sam.l - ref.l) * (sam.l - ref.l) +
		     (sam.a - ref.a) * (sam.a - ref.a) +
		     (sam.b - ref.b) * (sam.b - ref.b));
}

static inline double delta_cie76_diff_fast(const struct cieLAB sam,
					   const struct cieLAB ref) {
	return sqrt(fma(sam.l - ref.l, sam.l - ref.l,
			 fma(sam.a - ref.a, sam.a - ref.a,
			     fma(sam.b - ref.b, sam.b - ref.b, 0))));
}

static inline double delta_cie94_diff_fast(const struct cieLAB sam,
					   const struct cieLAB ref) {
	const double K1 = 0.045;
	const double K2 = 0.015;

	const double dl = sam.l - ref.l;
	const double da = sam.a - ref.a;
	const double db = sam.b - ref.b;

	const double C1 = sqrt(square(sam.a) + square(sam.b));
	const double C2 = sqrt(square(ref.a) + square(ref.b));
	const double dC = C1 - C2;

	const double dH2 = square(da) + square(db) - square(dC);

	const double sl = 1.0;
	const double sc = fma(K1, C1, 1.0);
	const double sh = fma(K2, C1, 1.0);

	const double term_L = dl / sl;
	const double term_C = dC / sc;
	const double term_H = sqrt(dH2) / sh;

	return sqrt(fma(term_L, term_L, fma(term_C, term_C, term_H * term_H)));
}

static inline double delta_ciede2000_diff_fast(const struct cieLAB sam,
					       const struct cieLAB ref) {
	double L1 = sam.l, a1 = sam.a, b1 = sam.b;
	double L2 = ref.l, a2 = ref.a, b2 = ref.b;

	double C1   = sqrt(fma(a1, a1, b1 * b1));
	double C2   = sqrt(fma(a2, a2, b2 * b2));
	double avgC = 0.5 * (C1 + C2);

	double G   = 0.5 * (1.0 - sqrt((avgC * avgC * avgC) /
					((avgC * avgC * avgC) + 15625.0)));
	double a1p = fma(a1, 1.0 + G, 0.0);
	double a2p = fma(a2, 1.0 + G, 0.0);

	double C1p   = sqrt(fma(a1p, a1p, b1 * b1));
	double C2p   = sqrt(fma(a2p, a2p, b2 * b2));
	double avgCp = 0.5 * (C1p + C2p);

	double h1p = atan2(b1, a1p);
	double h2p = atan2(b2, a2p);
	if (h1p < 0.0)
		h1p += 2.0 * M_PI;
	if (h2p < 0.0)
		h2p += 2.0 * M_PI;

	double deltLp = L2 - L1;
	double deltCp = C2p - C1p;

	double dhp;
	if (fabs(h2p - h1p) <= M_PI)
		dhp = h2p - h1p;
	else
		dhp = (h2p <= h1p) ? h2p - h1p + 2.0 * M_PI
				   : h2p - h1p - 2.0 * M_PI;

	double deltHp = 2.0 * sqrt(C1p * C2p) * sin(dhp / 2.0);

	double avgLp = 0.5 * (L1 + L2);
	double avgHp = (fabs(h1p - h2p) > M_PI)
			   ? fmod((h1p + h2p + 2.0 * M_PI), 2.0 * M_PI) * 0.5
			   : 0.5 * (h1p + h2p);

	double T = 1.0 - 0.17 * cos(avgHp - M_PI / 6.0) +
		   0.24 * cos(2.0 * avgHp) +
		   0.32 * cos(3.0 * avgHp + M_PI / 30.0) -
		   0.20 * cos(4.0 * avgHp - (7.0 * M_PI / 20.0));

	double angle	   = (avgHp * 180.0 / M_PI - 275.0) / 25.0;
	double delta_theta = (M_PI / 6.0) * exp(-(angle * angle));

	double Rc = 2.0 * sqrt((avgCp * avgCp * avgCp) /
				((avgCp * avgCp * avgCp) + 15625.0));
	double Sl = 1.0 + (0.015 * ((avgLp - 50.0) * (avgLp - 50.0))) /
			      sqrt(20.0 + ((avgLp - 50.0) * (avgLp - 50.0)));
	double Sc = fma(0.045, avgCp, 1.0);
	double Sh = fma(0.015, avgCp * T, 1.0);
	double Rt = -sin(2.0 * delta_theta) * Rc;

	double termL = deltLp / Sl;
	double termC = deltCp / Sc;
	double termH = deltHp / Sh;

	return sqrt(fma(termL, termL, fma(termC, termC, termH * termH)) +
		     Rt * termC * termH);
}

double delta_ok_diff(struct Color *sam, struct Color *ref) {
	if (!Color_has_space(*sam, COLOR_OKLAB)) {
		convert_srgb_to_oklab(sam);
	}
	if (!Color_has_space(*ref, COLOR_OKLAB)) {
		convert_srgb_to_oklab(ref);
	}
	return delta_ok_diff_fast(sam->oklab, ref->oklab);
}

double delta_cie76_diff(struct Color *sam, struct Color *ref) {
	if (!Color_has_space(*sam, COLOR_CIELAB)) {
		convert_srgb_to_cielab(sam);
	}
	if (!Color_has_space(*ref, COLOR_CIELAB)) {
		convert_srgb_to_cielab(ref);
	}
	return delta_cie76_diff_fast(sam->cielab, ref->cielab);
}

double delta_cie94_diff(struct Color *sam, struct Color *ref) {
	if (!Color_has_space(*sam, COLOR_CIELAB)) {
		convert_srgb_to_cielab(sam);
	}
	if (!Color_has_space(*ref, COLOR_CIELAB)) {
		convert_srgb_to_cielab(ref);
	}

	return delta_cie94_diff_fast(sam->cielab, ref->cielab);
}

double delta_ciede2000_diff(struct Color *sam, struct Color *ref) {
	if (!Color_has_space(*sam, COLOR_CIELAB)) {
		convert_srgb_to_cielab(sam);
	}
	if (!Color_has_space(*ref, COLOR_CIELAB)) {
		convert_srgb_to_cielab(ref);
	}
	return delta_ciede2000_diff_fast(sam->cielab, ref->cielab);
}

void Color_calc_spaces(struct Color *color) {
	convert_srgb_to_oklab(color);
	convert_srgb_to_cielab(color);
	convert_srgb_to_grayscale(color);
}
