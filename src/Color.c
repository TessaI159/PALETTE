#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "Color.h"

#ifdef PALETTE_DEBUG
#include "test_shared.h"
#endif

#define COLORSPACE_COUNT 3
#define SQUARE(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))
#define SECOND_SURSOLID(x) ((x) * (x) * (x) * (x) * (x) * (x) * (x))
#define DEG2RAD(x) (((x) / (180.0)) * (M_PI))
#define RAD2DEG(x) (((x) / M_PI) * (180.0))

/* Constants for a standard D65/2 illuminant */
static const double X2	  = 0.95047;
static const double Y2	  = 1.000;
static const double Z2	  = 1.08883;
static const double DELTA = (6.0 / 29.0) * (6.0 / 29.0) * (6.0 / 29.0);
#ifndef M_PI
static const double M_PI = 3.14159265358979323846;
extern E2000_diff   e2000_diffs_calc[NUM_E2000_PAIR]
#endif

    /* TODO (Tess): Allow creation of colors with oklab */
    /* TODO (Tess): oklab to srgb conversion */

    static inline float
    linearize(float channel) {
	return (channel > 0.04045f) ? powf((channel + 0.055f) / 1.055f, 2.4f)
				    : channel / 12.92f;
}

/* This is a special tool we will need for later */
static inline float gamma_correct(float channel) {
	return (channel > 0.0031308f)
		   ? 1.055f * powf(channel, 1.0f / 2.4f) - 0.055f
		   : 12.92f * channel;
}

struct Color Color_create(const float r, const float g, const float b) {
	struct Color color = {0};
	if (r < 0.0f || r > 255.0f || g < 0.0f || g > 255.0f || b < 0.0f ||
	    b > 255.0f) {
		fprintf(stderr, "rgb values out of range.\n");
	}

	color.data.srgb.r   = linearize(r / 255.0f);
	color.data.srgb.g   = linearize(g / 255.0f);
	color.data.srgb.b   = linearize(b / 255.0f);
	color.current_space = COLOR_SRGB;
	return color;
}

struct Color Color_create_norm(const float r, const float g, const float b) {
	struct Color color;
	color.data.srgb.r   = linearize(r);
	color.data.srgb.g   = linearize(g);
	color.data.srgb.b   = linearize(b);
	color.current_space = COLOR_SRGB;
	return color;
}

struct Color Color_create_lab(const float l, const float a, const float b) {
	struct Color color;
	color.data.cielab.l = l;
	color.data.cielab.a = a;
	color.data.cielab.b = b;
	color.current_space = COLOR_CIELAB;
	return color;
}

struct Color Color_create_ok(const float l, const float a, const float b) {
	struct Color color;
	color.data.oklab.l  = l;
	color.data.cielab.a = a;
	color.data.cielab.b = b;
	color.current_space = COLOR_CIELAB;
	return color;
}

static inline void convert_cielab_to_srgb(struct Color *color) {
	const struct cieLAB lab = color->data.cielab;

	const double fy = (lab.l + 16.0) / 116.0;
	const double fx = fy + (lab.a / 500.0);
	const double fz = fy - (lab.b / 200.0);

	double x = CUBE(fx) > DELTA ? CUBE(fx) : (fx - 16.0 / 116.0) / 7.787;
	double y =
	    lab.l > (DELTA * 116.0) ? CUBE(fy) : (fy - 16.0 / 116.0) / 7.787;
	double z = CUBE(fz) > DELTA ? CUBE(fz) : (fz - 16.0 / 116.0) / 7.787;

	x *= X2;
	y *= Y2;
	z *= Z2;

	color->data.srgb.r =
	    fma(x, 3.2404542, fma(y, -1.5371385, z * -0.4985314));
	color->data.srgb.g =
	    fma(x, -0.9692660, fma(y, 1.8760108, z * 0.0415560));
	color->data.srgb.b =
	    fma(x, 0.0556434, fma(y, -0.2040259, z * 1.0572252));

	color->current_space = COLOR_SRGB;
}

static inline void convert_oklab_to_srgb(struct Color *color) {
	const double L = color->data.oklab.l;
	const double A = color->data.oklab.a;
	const double B = color->data.oklab.b;

	const double l_ = fma(0.3963377774, A, fma(0.2158037573, B, L));
	const double m_ = fma(-0.1055613458, A, fma(-0.0638541728, B, L));
	const double s_ = fma(-0.0894841775, A, fma(-1.2914855480, B, L));

	const double l3 = l_ * l_ * l_;
	const double m3 = m_ * m_ * m_;
	const double s3 = s_ * s_ * s_;

	color->data.srgb.r =
	    fma(4.0767416621, l3, fma(-3.3077115913, m3, 0.2309699292 * s3));
	color->data.srgb.g =
	    fma(-1.2684380046, l3, fma(2.6097574011, m3, -0.3413193965 * s3));
	color->data.srgb.b =
	    fma(-0.0041960863, l3, fma(-0.7034186147, m3, 1.7076147010 * s3));

	color->current_space = COLOR_SRGB;
}

static void convert_srgb_to_cielab(struct Color *color) {
	const struct sRGB srgb = color->data.srgb;

	const double r = srgb.r;
	const double g = srgb.g;
	const double b = srgb.b;

	double x = fma(r, 0.4124564, fma(g, 0.3575761, b * 0.1804375));
	double y = fma(r, 0.2126729, fma(g, 0.7151522, b * 0.0721750));
	double z = fma(r, 0.0193339, fma(g, 0.1191920, b * 0.9503041));

	x /= X2;
	y /= Y2;
	z /= Z2;

	const double fx =
	    x > DELTA ? pow(x, 1.0 / 3.0) : fma(7.787, x, 16.0 / 116.0);
	const double fy =
	    y > DELTA ? pow(y, 1.0 / 3.0) : fma(7.787, y, 16.0 / 116.0);
	const double fz =
	    z > DELTA ? pow(z, 1.0 / 3.0) : fma(7.787, z, 16.0 / 116.0);

	color->data.cielab.l = fma(116.0, fy, -16.0);
	color->data.cielab.a = 500.0 * (fx - fy);
	color->data.cielab.b = 200.0 * (fy - fz);

	color->current_space = COLOR_CIELAB;
}

#ifdef PALETTE_DEBUG
void convert_srgb_to_oklab(struct Color *color) {
#else
static void convert_srgb_to_oklab(struct Color *color) {
#endif
	const struct sRGB srgb = color->data.srgb;

	double l = fma(0.4122214708, srgb.r,
		       fma(0.5363325363, srgb.g, 0.0514459929 * srgb.b));
	double m = fma(0.2119034982, srgb.r,
		       fma(0.6806995451, srgb.g, 0.1073969566 * srgb.b));
	double s = fma(0.0883024619, srgb.r,
		       fma(0.2817188376, srgb.g, 0.6299787005 * srgb.b));

	l = pow(l, 1.0 / 3.0);
	m = pow(m, 1.0 / 3.0);
	s = pow(s, 1.0 / 3.0);

	color->data.oklab.l =
	    fma(0.2104542553, l, fma(0.7936177850, m, -0.0040720468 * s));
	color->data.oklab.a =
	    fma(1.9779984951, l, fma(-2.4285922050, m, 0.4505937099 * s));
	color->data.oklab.b =
	    fma(0.0259040371, l, fma(0.7827717662, m, -0.8086757660 * s));

	color->current_space = COLOR_OKLAB;
}

static inline void convert_cielab_to_oklab(struct Color *color) {
	convert_cielab_to_srgb(color);
	convert_srgb_to_oklab(color);
}

static inline void convert_oklab_to_cielab(struct Color *color) {
	convert_oklab_to_srgb(color);
	convert_srgb_to_cielab(color);
}

typedef void (*convert_fn)(struct Color *);

convert_fn conversion_table[COLORSPACE_COUNT][COLORSPACE_COUNT] = {
    [COLOR_SRGB]   = {[COLOR_SRGB]   = NULL,
		      [COLOR_CIELAB] = convert_srgb_to_cielab,
		      [COLOR_OKLAB]  = convert_srgb_to_oklab},
    [COLOR_CIELAB] = {[COLOR_CIELAB] = NULL,
		      [COLOR_OKLAB]  = convert_cielab_to_oklab,
		      [COLOR_SRGB]   = convert_cielab_to_srgb},
    [COLOR_OKLAB]  = {[COLOR_OKLAB]  = NULL,
		      [COLOR_SRGB]   = convert_oklab_to_srgb,
		      [COLOR_CIELAB] = convert_oklab_to_cielab}};

#ifdef PALETTE_DEBUG
void convert_to(struct Color *color, enum ColorSpace t) {
#else
static inline void convert_to(struct Color *color, enum ColorSpace t) {
#endif

	if (t == color->current_space) {
		return;
	}

	convert_fn fn = conversion_table[color->current_space][t];
	if (fn) {
		fn(color);
	} else {
		printf("Unknown color conversion.\n");
	}
}

/* Color difference functions begin here */
static inline float delta_ok_diff_fast(const struct okLAB *sam,
				       const struct okLAB *ref) {
	return sqrtf(SQUARE(sam->l - ref->l) + SQUARE(sam->a - ref->a) +
		     SQUARE(sam->b - ref->b));
}

static inline float delta_cie76_diff_fast(const struct cieLAB *sam,
					  const struct cieLAB *ref) {
	return sqrtf(SQUARE(sam->l - ref->l) + SQUARE(sam->a - ref->a) +
		     SQUARE(sam->b - ref->b));
}

static inline float delta_cie94_diff_fast(const struct cieLAB *sam,
					  const struct cieLAB *ref) {
	struct cieLAB s = *sam;
	struct cieLAB r = *ref;

	const double K1 = 0.045;
	const double K2 = 0.015;

	const double dl = s.l - r.l;
	const double da = s.a - r.a;
	const double db = s.b - r.b;

	const double C1 = sqrt(SQUARE(s.a) + SQUARE(s.b));
	const double C2 = sqrt(SQUARE(r.a) + SQUARE(r.b));
	const double dC = C1 - C2;

	const double dH2 = SQUARE(da) + SQUARE(db) - SQUARE(dC);

	const double sl = 1.0;
	const double sc = fma(K1, C1, 1.0);
	const double sh = fma(K2, C1, 1.0);

	const double term_L = dl / sl;
	const double term_C = dC / sc;
	const double term_H = sqrt(dH2 < 0.0 ? 0.0 : dH2) / sh;

	return sqrt(SQUARE(term_L) + SQUARE(term_C) + SQUARE(term_H));
}

#ifdef PALETTE_DEBUG
float delta_ciede2000_diff_fast(E2000_diff *diff) {
#else
static inline float delta_ciede2000_diff_fast(const struct cieLAB *sam,
					      const struct cieLAB *ref) {
#endif

	const double kl = 1;
	const double kh = 1;
	const double kc = 1;
#ifdef PALETTE_DEBUG
	const double l1 = diff->l[0], a1 = diff->a[0], b1 = diff->b[0];
	const double l2 = diff->l[1], a2 = diff->a[1], b2 = diff->b[1];
#else
	const double l1 = sam->l, a1 = sam->a, b1 = sam->b;
	const double l2 = ref->l, a2 = ref->a, b2 = ref->b;
#endif

	/* Step 1 */
	const double c1 = sqrtf(SQUARE(a1) + SQUARE(b1));
	const double c2 = sqrtf(SQUARE(a2) + SQUARE(b2));

	const double avgc = (c1 + c2) / 2.0;
	/* Saves calculating this twice */
	const double g =
	    0.5f * (1.0f - sqrtf(SECOND_SURSOLID(avgc) /
				 (SECOND_SURSOLID(avgc) + 6103515625.0)));
	const double a1p = (1.0 + g) * a1;
	const double a2p = (1.0 + g) * a2;
	const double c1p = sqrt(SQUARE(a1p) + SQUARE(b1));
	const double c2p = sqrt(SQUARE(a2p) + SQUARE(b2));
	const double h1p = (b1 == 0.0 && a1p == 0.0)
			       ? 0
			       : fmod(RAD2DEG(atan2(b1, a1p)) + 360.0, 360.0);
	const double h2p = (b2 == 0.0f && a2p == 0.0f)
			       ? 0
			       : fmod(RAD2DEG(atan2(b2, a2p)) + 360.0, 360.0);

	/* Step 2 */
	const double dlp = l2 - l1;
	const double dcp = c2p - c1p;

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

	const double dHp = 2.0 * sqrt(c1p * c2p) * sin(DEG2RAD(dhp / 2.0));

	/* Step 3 */
	const double avglp = (l1 + l2) / 2.0;
	const double avgcp = (c1p + c2p) / 2.0;

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

	const double t = 1.0 - 0.17 * cos(DEG2RAD(avghp - 30.0)) +
			 0.24 * cos(DEG2RAD(2.0 * avghp)) +
			 0.32 * cos(DEG2RAD(3.0 * avghp + 6.0)) -
			 0.20 * cos(DEG2RAD(4.0 * avghp - 63.0));
	const double dtheta = 30.0 * exp(-(SQUARE((avghp - 275.0) / 25.0)));
	const double rc	    = 2.0 * sqrt(SECOND_SURSOLID(avgcp) /
					 (SECOND_SURSOLID(avgcp) + 6103515625.0));
	const double sl	    = 1.0 + ((0.015 * (SQUARE(avglp - 50.0))) /
				     sqrt(20.0 + (SQUARE(avglp - 50.0))));
	const double sc	    = 1.0 + 0.045 * avgcp;
	const double sh	    = 1.0 + 0.015 * avgcp * t;
	const double rt	    = -sin(DEG2RAD(2.0 * dtheta)) * rc;
	const double terml  = dlp / (kl * sl);
	const double termc  = dcp / (kc * sc);
	const double termh  = dHp / (kh * sh);

	const double ret = sqrt(SQUARE(terml) + SQUARE(termc) + SQUARE(termh) +
				rt * termc * termh);
#ifdef PALETTE_DEBUG
	diff->ap[0] = a1p;
	diff->ap[1] = a2p;
	diff->cp[0] = c1p;
	diff->cp[1] = c2p;
	diff->hp[0] = h1p;
	diff->hp[1] = h2p;

	diff->avghp = avghp;
	diff->g	    = g;
	diff->t	    = t;
	diff->sl    = sl;
	diff->sc    = sc;
	diff->sh    = sh;
	diff->rt    = rt;
	diff->diff  = ret;
#endif

	return ret;
}

float delta_ok_diff(struct Color *sam, struct Color *ref) {
	convert_to(sam, COLOR_OKLAB);
	convert_to(ref, COLOR_OKLAB);
	return delta_ok_diff_fast(&sam->data.oklab, &ref->data.oklab);
}

float delta_cie76_diff(struct Color *sam, struct Color *ref) {
	convert_to(sam, COLOR_CIELAB);
	convert_to(ref, COLOR_CIELAB);
	return delta_cie76_diff_fast(&sam->data.cielab, &ref->data.cielab);
}

float delta_cie94_diff(struct Color *sam, struct Color *ref) {
	convert_to(sam, COLOR_CIELAB);
	convert_to(ref, COLOR_CIELAB);
	return delta_cie94_diff_fast(&sam->data.cielab, &ref->data.cielab);
}

float delta_ciede2000_diff(struct Color *sam, struct Color *ref) {
	convert_to(sam, COLOR_CIELAB);
	convert_to(ref, COLOR_CIELAB);
#ifdef PALETTE_DEBUG
	struct E2000_diff diff;
	diff.l[0] = sam->data.cielab.l;
	diff.a[0] = sam->data.cielab.a;
	diff.b[0] = sam->data.cielab.b;
	diff.l[1] = ref->data.cielab.l;
	diff.a[1] = ref->data.cielab.a;
	diff.b[1] = ref->data.cielab.b;
	return delta_ciede2000_diff_fast(&diff);
#else
	return delta_ciede2000_diff_fast(&sam->data.cielab, &ref->data.cielab);
#endif
}

#ifdef PALETTE_DEBUG
void Color_print(struct Color *color) {
	convert_to(color, COLOR_SRGB);
	printf("Linear sRGB: (%f, %f, %f)\n", color->data.srgb.r,
	       color->data.srgb.g, color->data.srgb.b);
	convert_to(color, COLOR_CIELAB);
	printf("cieLAB: (%f, %f, %f)\n", color->data.cielab.l,
	       color->data.cielab.a, color->data.cielab.b);
	convert_to(color, COLOR_OKLAB);
	printf("okLAB: (%f, %f, %f)\n", color->data.oklab.l,
	       color->data.oklab.a, color->data.oklab.b);
}
#endif
