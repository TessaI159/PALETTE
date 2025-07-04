#define _USE_MATH_DEFINES
#include <math.h>

#include "Color.h"

#define COLOR_SPACE_BIT(space) (1u << (space))

/* Constants for a standard D65/2 illuminant */
static const float X2 = 95.047f;
static const float Y2 = 100.0f;
static const float Z2 = 108.883f;
#ifndef M_PI
static const float M_PI = 3.14159265358979323846f;
#endif

/* Mark/check the validity of colors */
static inline void Color_mark_space(struct Color   *color,
				    enum ColorSpace space) {
	color->valid_spaces |= COLOR_SPACE_BIT(space);
}

static inline int Color_has_space(const struct Color color,
				  enum ColorSpace    space) {
	return (color.valid_spaces & COLOR_SPACE_BIT(space)) != 0;
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

struct Color Color_create(const uint8_t r, const uint8_t g, const uint8_t b) {
	struct Color color;
	color.srgb.r = linearize((float)r);
	color.srgb.g = linearize((float)g);
	color.srgb.b = linearize((float)b);
	Color_mark_space(&color, COLOR_SRGB);
	return color;
}

static void convert_srgb_to_cielab(struct Color *color) {
	const struct sRGB srgb = color->srgb;

	float r = srgb.r * 100.0f;
	float g = srgb.g * 100.0f;
	float b = srgb.b * 100.0f;
	float x =
	    fmaf(r, 0.4124f, fmaf(g, 0.3576f, fmaf(b, 0.1805f, 0.0f))) / X2;
	float y =
	    fmaf(r, 0.2126f, fmaf(g, 0.7152f, fmaf(b, 0.0722f, 0.0f))) / Y2;
	float z =
	    fmaf(r, 0.0193f, fmaf(g, 0.1192f, fmaf(b, 0.9505f, 0.0f))) / Z2;
	double fx = (x > 0.008856f) ? cbrt(x) : fmaf(7.787f, x, 16.0f / 116.0f);
	double fy = (y > 0.008856f) ? cbrt(y) : fmaf(7.787f, y, 16.0f / 116.0f);
	double fz = (z > 0.008856f) ? cbrt(z) : fmaf(7.787f, z, 16.0f / 116.0f);

	struct cieLAB lab;
	lab.l = (y > 0.008856f) ? fmaf(116.0f, fy, -16.0f) : 903.3f * y;
	lab.a = fmaf(500.0f, fx - fy, 0.0f);
	lab.b = fmaf(200.0f, fy - fz, 0.0f);

	color->cielab = lab;
	Color_mark_space(color, COLOR_CIELAB);
}

static void convert_srgb_to_oklab(struct Color *color) {
	const struct sRGB srgb = color->srgb;

	float r = srgb.r * 100.0f;
	float g = srgb.g * 100.0f;
	float b = srgb.b * 100.0f;

	float x = fmaf(r, 0.4124f, fmaf(g, 0.3576f, fmaf(b, 0.1805f, 0.0f)));
	float y = fmaf(r, 0.2126f, fmaf(g, 0.7152f, fmaf(b, 0.0722f, 0.0f)));
	float z = fmaf(r, 0.0193f, fmaf(g, 0.1192f, fmaf(b, 0.9505f, 0.0f)));
	float l =
	    fmaf(0.8189330101f, x, fmaf(0.3618667424f, y, -0.1288597137f * z));
	float m =
	    fmaf(0.0329845436f, x, fmaf(0.9293118715f, y, 0.0361456387f * z));
	float s =
	    fmaf(0.0482003018f, x, fmaf(0.2643662691f, y, 0.6338517070f * z));

	l = cbrt(l);
	m = cbrt(m);
	s = cbrt(s);

	struct okLAB lab;
	lab.l = fmaf(0.2104542553f, l,
		     fmaf(0.7936177850f, m, fmaf(-0.0040720468f, s, 0.0f)));
	lab.a = fmaf(1.9779984951f, l,
		     fmaf(-2.4285922050f, m, fmaf(0.4505937099f, s, 0.0f)));
	lab.b = fmaf(0.0259040371f, l,
		     fmaf(0.7827717662f, m, fmaf(-0.8086757660f, s, 0.0f)));

	color->oklab = lab;
	Color_mark_space(color, COLOR_OKLAB);
}

static inline void convert_srgb_to_grayscale(struct Color *color) {
	struct sRGB srgb = color->srgb;

	float l = fmaf(0.2126f, srgb.r,
		       fmaf(0.7152f, srgb.g, fmaf(0.0722f, srgb.b, 0)));

	color->grayscale.l = l;
}

/* Color difference functions begin here */
static inline float euclidean_diff_fast(const struct sRGB sam,
					const struct sRGB ref) {
	return sqrtf(fmaf(sam.r - ref.r, sam.r - ref.r,
			  fmaf(sam.g - ref.g, sam.g - ref.g,
			       fmaf(sam.b - ref.b, sam.b - ref.b, 0))));
}

static inline float redmean_diff_fast(const struct sRGB sam,
				      const struct sRGB ref) {
	return fabs(sqrtf((2 + ((0.5f * (sam.r + ref.r)) / 256.0f)) *
			      (sam.r - ref.r) * (sam.r - ref.r) +
			  4 * ((sam.g - ref.g) * (sam.g - ref.g)) +
			  (2 + ((255.0f - (0.5f * (sam.r + ref.r))) / 256.0f)) *
			      (sam.b - ref.b) * (sam.b - ref.b)));
}

static inline float delta_ok_diff_fast(const struct okLAB sam,
				       const struct okLAB ref) {
	return sqrtf((sam.l - ref.l) * (sam.l - ref.l) +
		     (sam.a - ref.a) * (sam.a - ref.a) +
		     (sam.b - ref.b) * (sam.b - ref.b));
}

static inline float delta_cie76_diff_fast(const struct cieLAB sam,
					  const struct cieLAB ref) {
	return sqrtf(fmaf(sam.l - ref.l, sam.l - ref.l,
			  fmaf(sam.a - ref.a, sam.a - ref.a,
			       fmaf(sam.b - ref.b, sam.b - ref.b, 0))));
}

static inline float delta_cie94_diff_fast(const struct cieLAB sam,
					  const struct cieLAB ref) {
	const float K1 = 0.045f;
	const float K2 = 0.015f;

	const float dl = sam.l - ref.l;
	const float da = sam.a - ref.a;
	const float db = sam.b - ref.b;

	const float C1 = sqrtf(fmaf(sam.a, sam.a, sam.b * sam.b));
	const float C2 = sqrtf(fmaf(ref.a, ref.a, ref.b * ref.b));
	const float dC = C1 - C2;

	const float dH2 = fmaf(da, da, fmaf(db, db, 0.0f)) - fmaf(dC, dC, 0);

	const float sl = 1.0f;
	const float sc = fmaf(K1, C1, 1.0f);
	const float sh = fmaf(K2, C1, 1.0f);

	const float term_L = dl / sl;
	const float term_C = dC / sc;
	const float term_H = sqrtf(dH2) / sh;

	return sqrtf(
	    fmaf(term_L, term_L, fmaf(term_C, term_C, term_H * term_H)));
}

static inline float delta_ciede2000_diff_fast(const struct cieLAB sam,
					      const struct cieLAB ref) {
	float L1 = sam.l, a1 = sam.a, b1 = sam.b;
	float L2 = ref.l, a2 = ref.a, b2 = ref.b;

	float C1   = sqrtf(fmaf(a1, a1, b1 * b1));
	float C2   = sqrtf(fmaf(a2, a2, b2 * b2));
	float avgC = 0.5f * (C1 + C2);

	float G	  = 0.5f * (1.0f - sqrtf((avgC * avgC * avgC) /
					 ((avgC * avgC * avgC) + 15625.0f)));
	float a1p = fmaf(a1, 1.0f + G, 0.0f);
	float a2p = fmaf(a2, 1.0f + G, 0.0f);

	float C1p   = sqrtf(fmaf(a1p, a1p, b1 * b1));
	float C2p   = sqrtf(fmaf(a2p, a2p, b2 * b2));
	float avgCp = 0.5f * (C1p + C2p);

	float h1p = atan2f(b1, a1p);
	float h2p = atan2f(b2, a2p);
	if (h1p < 0.0f)
		h1p += 2.0f * M_PI;
	if (h2p < 0.0f)
		h2p += 2.0f * M_PI;

	float deltLp = L2 - L1;
	float deltCp = C2p - C1p;

	float dhp;
	if (fabs(h2p - h1p) <= M_PI)
		dhp = h2p - h1p;
	else
		dhp = (h2p <= h1p) ? h2p - h1p + 2.0f * M_PI
				   : h2p - h1p - 2.0f * M_PI;

	float deltHp = 2.0f * sqrtf(C1p * C2p) * sinf(dhp / 2.0f);

	float avgLp = 0.5f * (L1 + L2);
	float avgHp = (fabs(h1p - h2p) > M_PI)
			  ? fmod((h1p + h2p + 2.0f * M_PI), 2.0f * M_PI) * 0.5f
			  : 0.5f * (h1p + h2p);

	float T = 1.0f - 0.17f * cosf(avgHp - M_PI / 6.0f) +
		  0.24f * cosf(2.0f * avgHp) +
		  0.32f * cosf(3.0f * avgHp + M_PI / 30.0f) -
		  0.20f * cosf(4.0f * avgHp - (7.0f * M_PI / 20.0f));

	float delta_theta =
	    (M_PI / 6.0f) * exp(-(((avgHp * 180.0f / M_PI - 275.0f) / 25.0f) *
				  ((avgHp * 180.0f / M_PI - 275.0f / 25.0f))));
	float Rc = 2.0f * sqrtf((avgCp * avgCp * avgCp) /
				((avgCp * avgCp * avgCp) + 15625.0f));
	float Sl =
	    1.0f + (0.015f * ((avgLp - 50.0f) * (avgLp - 50.0f))) /
		       sqrtf(20.0f + ((avgLp - 50.0f) * (avgLp - 50.0f)));
	float Sc = fmaf(0.045f, avgCp, 1.0f);
	float Sh = fmaf(0.015f, avgCp * T, 1.0f);
	float Rt = -sinf(2.0f * delta_theta) * Rc;

	float termL = deltLp / Sl;
	float termC = deltCp / Sc;
	float termH = deltHp / Sh;

	return sqrtf(fmaf(termL, termL, fmaf(termC, termC, termH * termH)) +
		     Rt * termC * termH);
}

float euclidean_diff(struct Color *sam, struct Color *ref) {
	return euclidean_diff_fast(sam->srgb, ref->srgb);
}

float redmean_diff(struct Color *sam, struct Color *ref) {
	return redmean_diff_fast(sam->srgb, ref->srgb);
}

float delta_ok_diff(struct Color *sam, struct Color *ref) {
	if (!Color_has_space(*sam, COLOR_OKLAB)) {
		convert_srgb_to_oklab(sam);
	}
	if (!Color_has_space(*ref, COLOR_OKLAB)) {
		convert_srgb_to_oklab(ref);
	}
	return delta_ok_diff_fast(sam->oklab, ref->oklab);
}

float delta_cie76_diff(struct Color *sam, struct Color *ref) {
	if (!Color_has_space(*sam, COLOR_CIELAB)) {
		convert_srgb_to_cielab(sam);
	}
	if (!Color_has_space(*ref, COLOR_CIELAB)) {
		convert_srgb_to_cielab(ref);
	}
	return delta_cie76_diff_fast(sam->cielab, ref->cielab);
}

float delta_cie94_diff(struct Color *sam, struct Color *ref) {
	if (!Color_has_space(*sam, COLOR_CIELAB)) {
		convert_srgb_to_cielab(sam);
	}
	if (!Color_has_space(*ref, COLOR_CIELAB)) {
		convert_srgb_to_cielab(ref);
	}

	return delta_cie94_diff_fast(sam->cielab, ref->cielab);
}

float delta_ciede2000_diff(struct Color *sam, struct Color *ref) {
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
