/*
 * Thanks to easyrgb.com Björn Ottosson and Wikipedia for the color conversion
 * and difference formulas
 * https://www.easyrgb.com/en/math.php
 * en.wikipedia.org/wiki/Color_difference
 * en.wikipedia.org/wiki/Oklab_color_space
 * bottosson.github.io/posts/oklab
 */

#include "Color.h"
#include <math.h>
#include <stdalign.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define COLOR_SPACE_BIT(space) (1u << (space))
#define square(num) ((num) * (num))
#define cube(num) ((num) * (num) * (num))

/*
 * Constants for a standard D65/2 illuminant
 * https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D
 * https://en.wikipedia.org/wiki/CIE_1931_color_space#CIE_standard_observer
 */

static const double X2 = 95.047;
static const double Y2 = 100.0;
static const double Z2 = 108.883;
#ifndef M_PI
static const double M_PI = 3.14159265358979323846;
#endif
/* static const double EPSILON = 1e-5; */

enum ColorSpace {
	COLOR_OKLAB = 1,
	COLOR_CIELAB,
	COLOR_SRGB,
	COLOR_XYZ,
};

/* Supported color spaces */
struct okLAB {
	double l, a, b;
};

struct sRGB {
	double r, g, b;
};

struct XYZ {
	double x, y, z;
};

struct cieLAB {
	double l, a, b;
};

struct Color {
	struct okLAB oklab;
	struct cieLAB cielab;
	struct sRGB srgb;
	struct XYZ xyz;
	uint32_t valid_spaces;
} __attribute__((aligned(64)));

uint64_t Color_size() {
	return sizeof(Color);
}

Color* Color_create(uint8_t r, uint8_t g, uint8_t b) {
	Color* rt = malloc(sizeof(Color));
	if (!rt) {
		perror("Color malloc failed.\n");
		return NULL;
	}
	memset(rt, 0, sizeof(Color));
	rt->srgb.r = r;
	rt->srgb.g = g;
	rt->srgb.b = b;

	return rt;
}

void Color_destroy(Color* color) {
	if (!color)
		return;
	free(color);
}

/* Mark/check the validity of colors */
static inline void Color_mark_space(Color* color, enum ColorSpace space) {
	color->valid_spaces |= COLOR_SPACE_BIT(space);
}

static inline int Color_has_space(Color* color, enum ColorSpace space) {
	return (color->valid_spaces & COLOR_SPACE_BIT(space)) != 0;
}

/* These next functions are just helper functions for the color conversions */
static inline double linearize(double channel) {
	return (channel > 0.04045) ? pow((channel + 0.055) / 1.055, 2.4)
				   : channel / 12.92;
}

/* static inline double gamma_encode(double channel) { */
/* 	return channel <= 0.0031308 ? (12.92 * channel) */
/* 				    : (1.055 * pow(channel, 1.0 / 2.4) - 0.055);
 */
/* } */

/* static inline int almost_equal_double(double a, double b) { */
/* 	return fabs(a - b) < EPSILON; */
/* } */

/* static inline double max_double(double a, double b) { */
/* 	return a > b ? a : b; */
/* } */

/* static inline double min_double(double a, double b) { */
/* 	return a < b ? a : b; */
/* } */

/* static inline double clamp_double(double arg, double low, double high) { */
/* 	if (almost_equal_double(low, high)) { */
/* 		return low; */
/* 	} else if (low > high) { */
/* 		double temp = low; */
/* 		low	    = high; */
/* 		high	    = temp; */
/* 	} */
/* 	if (arg < low) { */
/* 		return low; */
/* 	} else if (arg > high) { */
/* 		return high; */
/* 	} */
/* 	return arg; */
/* } */

/* static inline double hue_to_rgb(double p, double q, double t) { */
/* 	if (t < 0.0) */
/* 		t += 1.0; */
/* 	if (t > 1.0) */
/* 		t -= 1.0; */
/* 	if (t < 1.0 / 6.0) */
/* 		return p + (q - p) * 6.0 * t; */
/* 	if (t < 1.0 / 2.0) */
/* 		return q; */
/* 	if (t < 2.0 / 3.0) */
/* 		return fma(p, q - p, (2.0 / 3.0 - t) * 6.0); */
/* 	return p; */
/* } */

/* Color conversion functions begin here */

static void convert_srgb_to_xyz(Color* color) {
	/* srgb is always valid, as colors as required to be created with srgb
	 * values */
	const struct sRGB srgb = color->srgb;

	double r = linearize(srgb.r / 255.0) * 100.0;
	double g = linearize(srgb.g / 255.0) * 100.0;
	double b = linearize(srgb.b / 255.0) * 100.0;

	struct XYZ xyz;
	xyz.x = fma(r, 0.4124, fma(g, 0.3576, fma(b, 0.1805, 0.0)));
	xyz.y = fma(r, 0.2126, fma(g, 0.7152, fma(b, 0.0722, 0.0)));
	xyz.z = fma(r, 0.0193, fma(g, 0.1192, fma(b, 0.9505, 0.0)));

	color->xyz = xyz;
	Color_mark_space(color, COLOR_XYZ);
}

/* static void convert_xyz_to_srgb(Color *color) { */
/* 	const struct XYZ xyz = color->xyz; */
/* 	double		 x   = xyz.x / 100.0; */
/* 	double		 y   = xyz.y / 100.0; */
/* 	double		 z   = xyz.z / 100.0; */
/* 	double r_lin = fma(x, 3.2406, fma(y, -1.5372, fma(z, -0.4986, 0.0))); */
/* 	double g_lin = fma(x, -0.9689, fma(y, 1.8758, fma(z, 0.0415, 0.0))); */
/* 	double b_lin = fma(x, 0.0557, fma(y, -0.2040, fma(z, 1.0570, 0.0))); */

/* 	struct sRGB srgb; */
/* 	srgb.r = round(clamp_double(gamma_encode(r_lin), 0.0, 1.0) * 255.0); */
/* 	srgb.g = round(clamp_double(gamma_encode(g_lin), 0.0, 1.0) * 255.0); */
/* 	srgb.b = round(clamp_double(gamma_encode(b_lin), 0.0, 1.0) * 255.0); */

/* 	color->srgb = srgb; */
/* 	Color_mark_space(color, COLOR_SRGB); */
/* } */

static void convert_xyz_to_cielab(Color* color) {
	const struct XYZ xyz = color->xyz;

	double x = xyz.x / X2;
	double y = xyz.y / Y2;
	double z = xyz.z / Z2;

	double fx = (x > 0.008856) ? cbrt(x) : fma(7.787, x, 16.0 / 116.0);
	double fy = (y > 0.008856) ? cbrt(y) : fma(7.787, y, 16.0 / 116.0);
	double fz = (z > 0.008856) ? cbrt(z) : fma(7.787, z, 16.0 / 116.0);

	struct cieLAB lab;
	lab.l = (y > 0.008856) ? fma(116.0, fy, -16.0) : 903.3 * y;
	lab.a = fma(500.0, fx - fy, 0.0);
	lab.b = fma(200.0, fy - fz, 0.0);

	color->cielab = lab;
	Color_mark_space(color,
			 COLOR_OKLAB); // Using okLAB slot to store CIELAB
}

/* static void convert_cielab_to_xyz(Color *color) { */
/* 	const struct cieLAB lab = color->cielab; */

/* 	double fy = fma(1.0 / 116.0, lab.l, 16.0 / 116.0); */
/* 	double fx = fma(1.0 / 500.0, lab.a, fy); */
/* 	double fz = fy - lab.b / 200.0; */

/* 	double fx3 = cube(fx); */
/* 	double fy3 = cube(fy); */
/* 	double fz3 = cube(fz); */

/* 	double xr = (fx3 > 0.008856) ? fx3 : (fx - 16.0 / 116.0) / 7.787; */
/* 	double yr = (lab.l > 7.9996) ? fy3 : lab.l / 903.3; */
/* 	double zr = (fz3 > 0.008856) ? fz3 : (fz - 16.0 / 116.0) / 7.787; */

/* 	struct XYZ xyz; */
/* 	xyz.x = X2 * xr; */
/* 	xyz.y = Y2 * yr; */
/* 	xyz.z = Z2 * zr; */

/* 	color->xyz = xyz; */
/* 	Color_mark_space(color, COLOR_XYZ); */
/* } */

static void convert_xyz_to_oklab(Color* color) {
	const struct XYZ xyz = color->xyz;
	double l = fma(0.8189330101, xyz.x,
		       fma(0.3618667424, xyz.y, -0.1288597137 * xyz.z));
	double m = fma(0.0329845436, xyz.x,
		       fma(0.9293118715, xyz.y, 0.0361456387 * xyz.z));
	double s = fma(0.0482003018, xyz.x,
		       fma(0.2643662691, xyz.y, 0.6338517070 * xyz.z));

	l = cbrt(l);
	m = cbrt(m);
	s = cbrt(s);

	struct okLAB lab;
	lab.l = fma(0.2104542553, l,
		    fma(0.7936177850, m, fma(-0.0040720468, s, 0.0)));
	lab.a = fma(1.9779984951, l,
		    fma(-2.4285922050, m, fma(0.4505937099, s, 0.0)));
	lab.b = fma(0.0259040371, l,
		    fma(0.7827717662, m, fma(-0.8086757660, s, 0.0)));

	color->oklab = lab;
	Color_mark_space(color, COLOR_OKLAB);
}

/* static void convert_oklab_to_xyz(Color *color) { */
/* 	const struct okLAB lab = color->oklab; */

/* 	double l = lab.l + 0.3963377774 * lab.a + 0.2158037573 * lab.b; */
/* 	double m = lab.l - 0.1055613458 * lab.a - 0.0638541728 * lab.b; */
/* 	double s = lab.l - 0.0894841775 * lab.a - 1.2914855480 * lab.b; */

/* 	l *= l * l; */
/* 	m *= m * m; */
/* 	s *= s * s; */

/* 	struct XYZ xyz; */
/* 	xyz.x	   = fma(1.2270138511, l, */
/* 			 fma(-0.5577999807, m, fma(0.2812561490, s, 0.0))); */
/* 	xyz.y	   = fma(-0.0405801784, l, */
/* 			 fma(1.1122568696, m, fma(-0.0716766787, s, 0.0))); */
/* 	xyz.z	   = fma(-0.0763812845, l, */
/* 			 fma(-0.4214819784, m, fma(1.5861632204, s, 0.0))); */
/* 	color->xyz = xyz; */
/* 	Color_mark_space(color, COLOR_XYZ); */
/* } */

/* static void convert_srgb_to_hsl(Color *color) { */
/* 	const struct sRGB srgb = color->srgb; */
/* 	double		  r    = srgb.r / 255.0; */
/* 	double		  g    = srgb.g / 255.0; */
/* 	double		  b    = srgb.b / 255.0; */

/* 	double max   = max_double(r, max_double(g, b)); */
/* 	double min   = min_double(r, min_double(g, b)); */
/* 	double delta = max - min; */

/* 	double h = 0.0; */
/* 	double s = 0.0; */
/* 	double l = (max + min) / 2.0; */

/* 	if (!almost_equal_double(delta, 0.0)) { */
/* 		if (almost_equal_double(max, r)) { */
/* 			h = 60.0 * (fmod((g - b) / delta, 6.0)); */
/* 		} else if (almost_equal_double(max, g)) { */
/* 			h = 60.0 * ((b - r) / delta + 2.0); */
/* 		} else { */
/* 			h = 60.0 * ((r - g) / delta + 4.0); */
/* 		} */

/* 		if (h < 0.0) */
/* 			h += 360.f; */

/* 		s = (l <= 0.5) ? delta / (max + min) */
/* 			       : delta / (2.0 - max - min); */
/* 	} */
/* 	struct HSL hsl; */
/* 	hsl.h = h; */
/* 	hsl.s = s * 100.0; */
/* 	hsl.l = l * 100.0; */

/* 	color->hsl = hsl; */
/* 	Color_mark_space(color, COLOR_HSL); */
/* } */

/* static void convert_hsl_to_srgb(Color *color) { */
/* 	const struct HSL hsl = color->hsl; */
/* 	double		 h   = hsl.h; */
/* 	double		 s   = hsl.s / 100.0; */
/* 	double		 l   = hsl.l / 100.0; */

/* 	double r; */
/* 	double g; */
/* 	double b; */

/* 	if (almost_equal_double(s, 0.0)) { */
/* 		r = g = b = l; */
/* 	} else { */
/* 		double q   = (l < 0.5) ? (l * (1.0 + s)) : (l + s - l * s); */
/* 		double p   = fma(2.0, l, -q); */
/* 		double h_k = h / 360.0; */
/* 		r	   = hue_to_rgb(p, q, h_k + 1.0 / 3.0); */
/* 		g	   = hue_to_rgb(p, q, h_k); */
/* 		b	   = hue_to_rgb(p, q, h_k - 1.0 / 3.0); */
/* 	} */

/* 	struct sRGB rgb; */
/* 	rgb.r	    = round(clamp_double(r, 0.0, 1.0) * 255.0); */
/* 	rgb.g	    = round(clamp_double(g, 0.0, 1.0) * 255.0); */
/* 	rgb.b	    = round(clamp_double(b, 0.0, 1.0) * 255.0); */
/* 	color->srgb = rgb; */
/* 	Color_mark_space(color, COLOR_SRGB); */
/* } */

/* static void convert_srgb_to_hsv(Color *color) { */
/* 	const struct sRGB srgb = color->srgb; */
/* 	double		  r    = srgb.r / 255.0; */
/* 	double		  g    = srgb.g / 255.0; */
/* 	double		  b    = srgb.b / 255.0; */

/* 	double max   = max_double(r, max_double(g, b)); */
/* 	double min   = min_double(r, min_double(g, b)); */
/* 	double delta = max - min; */

/* 	struct HSV hsv; */
/* 	hsv.v = max * 100.0; */
/* 	hsv.s = (almost_equal_double(max, 0.0)) ? 0.0 : (delta / max) * 100.0;
 */

/* 	if (almost_equal_double(delta, 0.0)) { */
/* 		hsv.h = 0.0; */
/* 	} else if (almost_equal_double(max, r)) { */
/* 		hsv.h = 60.0 * fmod((g - b) / delta, 6.0); */
/* 	} else if (almost_equal_double(max, g)) { */
/* 		hsv.h = 60.0 * ((b - r) / delta + 2.0); */
/* 	} else { */
/* 		hsv.h = 60.0 * ((r - g) / delta + 4.0); */
/* 	} */

/* 	if (hsv.h < 0.0) { */
/* 		hsv.h += 360; */
/* 	} */

/* 	color->hsv = hsv; */
/* 	Color_mark_space(color, COLOR_HSV); */
/* } */

/* static void convert_hsv_to_srgb(Color *color) { */
/* 	const struct HSV hsv = color->hsv; */
/* 	double		 h   = hsv.h; */
/* 	double		 s   = hsv.s / 100.0; */
/* 	double		 v   = hsv.v / 100.0; */

/* 	double c = v * s; */
/* 	double x = c * (1.0 - fabs(fmod(h / 60.0, 2) - 1)); */
/* 	double m = v - c; */

/* 	double r; */
/* 	double g; */
/* 	double b; */

/* 	if (h < 60.0) { */
/* 		r = c; */
/* 		g = x; */
/* 		b = 0; */
/* 	} else if (h < 120.0) { */
/* 		r = x; */
/* 		g = c; */
/* 		b = 0; */
/* 	} else if (h < 180.0) { */
/* 		r = 0; */
/* 		g = c; */
/* 		b = x; */
/* 	} else if (h < 240.0) { */
/* 		r = 0; */
/* 		g = x; */
/* 		b = c; */
/* 	} else if (h < 300.0) { */
/* 		r = x; */
/* 		g = 0; */
/* 		b = c; */
/* 	} else { */
/* 		r = c; */
/* 		g = 0; */
/* 		b = x; */
/* 	} */

/* 	struct sRGB rgb; */
/* 	rgb.r = round((r + m) * 255.0); */
/* 	rgb.g = round((g + m) * 255.0); */
/* 	rgb.b = round((b + m) * 255.0); */

/* 	color->srgb = rgb; */
/* 	Color_mark_space(color, COLOR_SRGB); */
/* } */

/* I know this is super fragile, but it's not meant to be extensible. Just
   useful */
static void convert_srgb_to_cielab(Color* color) {
	if (!Color_has_space(color, COLOR_XYZ)) {
		convert_srgb_to_xyz(color);
	}
	if (!Color_has_space(color, COLOR_CIELAB)) {
		convert_xyz_to_cielab(color);
	}
}

static void convert_srgb_to_oklab(Color* color) {
	if (!Color_has_space(color, COLOR_XYZ)) {
		convert_srgb_to_xyz(color);
	}
	if (!Color_has_space(color, COLOR_OKLAB)) {
		convert_xyz_to_oklab(color);
	}
}

/* Color difference functions begin here */

double euclidean_diff(Color* sam, Color* ref) {
	return sqrt(square(sam->srgb.r - ref->srgb.r) +
		    square(sam->srgb.g - ref->srgb.g) +
		    square(sam->srgb.b - ref->srgb.b));
}

double redmean_diff(Color* sam, Color* ref) {
	return fabs(
	    sqrt((2 + ((0.5 * (sam->srgb.r + ref->srgb.r)) / 256.0)) *
		     square((sam->srgb.r - ref->srgb.r)) +
		 4 * square((sam->srgb.g - ref->srgb.g)) +
		 (2 + ((255.0 - (0.5 * (sam->srgb.r + ref->srgb.r))) / 256.0)) *
		     square((sam->srgb.b - ref->srgb.b))));
}

double delta_ok_diff(Color* sam, Color* ref) {
	if (!Color_has_space(sam, COLOR_OKLAB)) {
		convert_srgb_to_oklab(sam);
	}
	if (!Color_has_space(ref, COLOR_OKLAB)) {
		convert_srgb_to_oklab(ref);
	}

	return sqrt(square(sam->oklab.l - ref->oklab.l) +
		    square(sam->oklab.a - ref->oklab.a) +
		    square(sam->oklab.b - ref->oklab.b));
}

double delta_cie76_diff(Color* sam, Color* ref) {
	if (!Color_has_space(sam, COLOR_CIELAB)) {
		convert_srgb_to_cielab(sam);
	}
	if (!Color_has_space(ref, COLOR_CIELAB)) {
		convert_srgb_to_cielab(ref);
	}
	return sqrt(fma(
	    sam->cielab.l - ref->cielab.l, sam->cielab.l - ref->cielab.l,
	    fma(sam->cielab.a - ref->cielab.a, sam->cielab.a - ref->cielab.a,
		sam->cielab.b - ref->cielab.b * sam->cielab.b -
		    ref->cielab.b)));
}

double delta_cie94_diff(Color* sam, Color* ref) {
	if (!Color_has_space(sam, COLOR_CIELAB)) {
		convert_srgb_to_cielab(sam);
	}
	if (!Color_has_space(ref, COLOR_CIELAB)) {
		convert_srgb_to_cielab(ref);
	}
	const double K1 = 0.045;
	const double K2 = 0.015;

	double dl = sam->cielab.l - ref->cielab.l;
	double da = sam->cielab.a - ref->cielab.a;
	double db = sam->cielab.b - ref->cielab.b;

	double C1 = sqrt(
	    fma(sam->cielab.a, sam->cielab.a, sam->cielab.b * sam->cielab.b));
	double C2 = sqrt(
	    fma(ref->cielab.a, ref->cielab.a, ref->cielab.b * ref->cielab.b));
	double dC = C1 - C2;

	double dH2 = fma(da, da, db * db) - dC * dC;

	double sl = 1.0;
	double sc = fma(K1, C1, 1.0);
	double sh = fma(K2, C1, 1.0);

	double term_L = dl / sl;
	double term_C = dC / sc;
	double term_H = sqrt(dH2) / sh;

	return sqrt(fma(term_L, term_L, fma(term_C, term_C, term_H * term_H)));
}

double delta_ciede2000_diff(Color* sam, Color* ref) {
	if (!Color_has_space(sam, COLOR_CIELAB)) {
		convert_srgb_to_cielab(sam);
	}
	if (!Color_has_space(ref, COLOR_CIELAB)) {
		convert_srgb_to_cielab(ref);
	}
	double L1 = sam->cielab.l, a1 = sam->cielab.a, b1 = sam->cielab.b;
	double L2 = ref->cielab.l, a2 = ref->cielab.a, b2 = ref->cielab.b;

	double C1 = sqrt(fma(a1, a1, b1 * b1));
	double C2 = sqrt(fma(a2, a2, b2 * b2));
	double avgC = 0.5 * (C1 + C2);

	double G = 0.5 * (1.0 - sqrt(cube(avgC) / (cube(avgC) + cube(25.0))));
	double a1p = fma(a1, 1.0 + G, 0.0);
	double a2p = fma(a2, 1.0 + G, 0.0);

	double C1p = sqrt(fma(a1p, a1p, b1 * b1));
	double C2p = sqrt(fma(a2p, a2p, b2 * b2));
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

	double delta_theta =
	    (M_PI / 6.0) * exp(-square((avgHp * 180.0 / M_PI - 275.0) / 25.0));
	double Rc = 2.0 * sqrt(cube(avgCp) / (cube(avgCp) + cube(25.0)));
	double Sl = 1.0 + (0.015 * square(avgLp - 50.0)) /
			      sqrt(20.0 + square(avgLp - 50.0));
	double Sc = fma(0.045, avgCp, 1.0);
	double Sh = fma(0.015, avgCp * T, 1.0);
	double Rt = -sin(2.0 * delta_theta) * Rc;

	double termL = deltLp / Sl;
	double termC = deltCp / Sc;
	double termH = deltHp / Sh;

	return sqrt(fma(termL, termL, fma(termC, termC, termH * termH)) +
		    Rt * termC * termH);
}

void Color_calc_spaces(Color* color) {
	convert_srgb_to_oklab(color);
	convert_srgb_to_cielab(color);
}
