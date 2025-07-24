#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#if defined(__x86_64__) || defined(_M_X64)
#include <pmmintrin.h>
#include <smmintrin.h>
#endif

#include "Color.h"

#define SQUARE_SCALAR(x) ((x) * (x))
#define POW_7_SCALAR(x) ((x) * (x) * (x) * (x) * (x) * (x) * (x))
#define DEG2RAD(x) (((x) / (180.0)) * (M_PI))
#define RAD2DEG(x) (((x) / M_PI) * (180.0))
#define _mm_abs_ps(x) _mm_and_ps((x), _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF)))
#define FMADD(a, b, c) _mm_add_ps(_mm_mul_ps((a), (b)), (c))
#define FL_PER_VEC 4
#define _MM_PI 3.14159265f
#define _MM_2PI 6.28318531f
#define _MM_PI_2 1.57079633f
#define _MM_1_OVER_2PI 0.159154943f

static inline __m128 pow7(__m128 x) {
	__m128 x2 = _mm_mul_ps(x, x);
	__m128 x4 = _mm_mul_ps(x2, x2);
	return _mm_mul_ps(x, _mm_mul_ps(x2, x4));
}

static inline __m128 square(__m128 x) {
	return _mm_mul_ps(x, x);
}

static inline __m128 deg_2_rad(__m128 x) {
	return _mm_mul_ps(_mm_div_ps(x, _mm_set1_ps(180.0f)),
			  _mm_set1_ps(_MM_PI));
}

static inline __m128 rad_2_deg(__m128 x) {
	return _mm_mul_ps(_mm_div_ps(x, _mm_set1_ps(_MM_PI)),
			  _mm_set1_ps(180.0f));
}

static inline __m128 V(float x) {
	return _mm_set1_ps(x);
}

static inline void
compute_tail(const struct Color *restrict colors,
	     const struct Color *restrict cents, float *diffs,
	     float (*scalar_diff)(const struct Color *restrict,
				  const struct Color *restrict, int, int)) {
	const uint32_t num_col	     = colors->num;
	const uint32_t num_cent	     = cents->num;
	size_t	       col_not_comp  = num_col - (num_col % FL_PER_VEC);
	size_t	       cent_not_comp = num_cent - (num_cent % FL_PER_VEC);
	for (size_t i = col_not_comp; i < num_col; ++i) {
		for (size_t j = 0; j < num_cent; ++j) {
			diffs[j * num_col + i] =
			    scalar_diff(colors, cents, i, j);
		}
	}
	for (size_t j = cent_not_comp; j < num_cent; ++j) {
		for (size_t i = 0; i < col_not_comp; ++i) {
			diffs[j * num_col + i] =
			    scalar_diff(colors, cents, i, j);
		}
	}
}

static inline __m128 _mm_atan_ps(__m128 x) {
	const __m128 one  = V(1.0f);
	const __m128 pi_2 = V(_MM_PI_2);

	__m128 sign_mask = _mm_and_ps(x, V(-0.0f));
	__m128 abs_x	 = _mm_andnot_ps(sign_mask, x);

	__m128 use_identity_mask = _mm_cmpgt_ps(abs_x, one);

	__m128 x_inv = _mm_div_ps(one, abs_x);

	__m128 x_red = _mm_or_ps(_mm_and_ps(use_identity_mask, x_inv),
				 _mm_andnot_ps(use_identity_mask, abs_x));

	__m128 x2 = _mm_mul_ps(x_red, x_red);
	__m128 y  = V(-0.01891985f);
	y	  = FMADD(y, x2, V(0.06909289f));
	y	  = FMADD(y, x2, V(-0.12945813f));
	y	  = FMADD(y, x2, V(0.19787976f));
	y	  = FMADD(y, x2, V(-0.33319414f));
	y	  = FMADD(y, x2, V(0.99999750f));

	y = _mm_mul_ps(y, x_red);

	__m128 atan_corrected = _mm_sub_ps(pi_2, y);
	y = _mm_or_ps(_mm_and_ps(use_identity_mask, atan_corrected),
		      _mm_andnot_ps(use_identity_mask, y));

	y = _mm_or_ps(y, sign_mask);
	return y;
}

static inline __m128 _mm_atan2_ps(__m128 y, __m128 x) {
	const __m128 zero     = V(0.0f);
	const __m128 pi	      = V(_MM_PI);
	const __m128 pi_2     = V(_MM_PI_2);
	const __m128 neg_pi_2 = V(-_MM_PI_2);

	__m128 abs_x = _mm_andnot_ps(V(-0.0f), x);

	__m128 safe_x = _mm_blendv_ps(x, V(1.0f), _mm_cmpeq_ps(abs_x, zero));
	__m128 z      = _mm_div_ps(y, safe_x);
	__m128 atan_z = _mm_atan_ps(z);

	__m128 x_lt0 = _mm_cmplt_ps(x, zero);
	__m128 y_eq0 = _mm_cmpeq_ps(y, zero);
	__m128 y_lt0 = _mm_cmplt_ps(y, zero);
	__m128 y_gt0 = _mm_cmpgt_ps(y, zero);
	__m128 x_eq0 = _mm_cmpeq_ps(x, zero);

	__m128 add_pi = _mm_and_ps(x_lt0, y_gt0);
	__m128 sub_pi = _mm_and_ps(x_lt0, y_lt0);
	__m128 result = atan_z;
	result	      = _mm_add_ps(result, _mm_and_ps(add_pi, pi));
	result	      = _mm_sub_ps(result, _mm_and_ps(sub_pi, pi));

	__m128 x0_y0   = _mm_and_ps(x_eq0, y_eq0);
	__m128 x0_ygt0 = _mm_and_ps(x_eq0, y_gt0);
	__m128 x0_ylt0 = _mm_and_ps(x_eq0, y_lt0);

	result = _mm_blendv_ps(result, zero, x0_y0);
	result = _mm_blendv_ps(result, pi_2, x0_ygt0);
	result = _mm_blendv_ps(result, neg_pi_2, x0_ylt0);

	return result;
}

static inline __m128 _mm_sin_ps(__m128 x) {
	const __m128 inv_2pi = V(_MM_1_OVER_2PI);
	const __m128 two_pi  = V(_MM_2PI);
	const __m128 pi	     = V(_MM_PI);

	__m128 k = _mm_mul_ps(x, inv_2pi);
	k	 = _mm_round_ps(k, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
	x	 = _mm_sub_ps(x, _mm_mul_ps(k, two_pi));
	x	 = _mm_sub_ps(x, pi);

	__m128 x2 = _mm_mul_ps(x, x);
	__m128 y  = V(0.00000227f);
	y	  = FMADD(y, x2, V(-0.00019506f));
	y	  = FMADD(y, x2, V(0.00832404f));
	y	  = FMADD(y, x2, V(-0.16665670f));
	y	  = FMADD(y, x2, V(0.99999708f));

	y = _mm_mul_ps(y, -x);
	return _mm_min_ps(_mm_max_ps(y, V(-1.0f)), V(1.0f));
}

static inline __m128 _mm_cos_ps(__m128 x) {
	const __m128 inv_2pi = V(_MM_1_OVER_2PI);
	const __m128 two_pi  = V(_MM_2PI);
	const __m128 pi	     = V(_MM_PI);

	__m128 k = _mm_mul_ps(x, inv_2pi);
	k	 = _mm_round_ps(k, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
	x	 = _mm_sub_ps(x, _mm_mul_ps(k, two_pi));
	x	 = _mm_sub_ps(x, pi);

	__m128 x2 = _mm_mul_ps(x, x);
	__m128 y  = V(0.00001974f);
	y	  = FMADD(y, x2, V(-0.00135731f));
	y	  = FMADD(y, x2, V(0.04159290f));
	y	  = FMADD(y, x2, V(-0.49994522f));
	y	  = FMADD(y, x2, V(0.99999571f));

	y = _mm_mul_ps(y, V(-1.0f));
	return _mm_min_ps(_mm_max_ps(y, V(-1.0f)), V(1.0f));
}

/* A note:
   This function is only tested on the interval [-87.33654, 88.72283], hence the
   clamping.
   But what if you pass a value outside that range, you say? Behold.

   30.0 * exp(-(SQUARE_SCALAR((avghp - 275.0) / 25.0)));
   This is the only line where we call exp in the entire translation unit.
   0 < avghp < 360, as avghp is an angle clamped to one rotation around the unit
   circle in degrees (ew).
   Then -11 <= (avghp - 275) / 25 <= 3.4 => exp(x) with -121 <= x <= 0
   |exp(-121) - exp(-87.33654)| = 1.1754e-38, which is well outside the bounds
   of either a float or a double (or even machine epsilon, for that matter).
   Even scaling it by a factor of thirty gives an error so miniscule it's
   unnoticeable.
 */
static inline __m128 _mm_exp_ps(__m128 x) {
	x	  = _mm_min_ps(_mm_max_ps(x, V(-87.33654)), V(88.72283));
	__m128	a = V(12102203.0f); /* (1 << 23) / log(2) */
	__m128i b = _mm_set1_epi32(127 * (1 << 23) - 298765);
	__m128i t = _mm_add_epi32(_mm_cvtps_epi32(_mm_mul_ps(a, x)), b);
	return _mm_castsi128_ps(t);
}

/* row * num_col + col */
/* Colors are columns, Centroids are rows */
static inline float delta_lab_euc_diff_scalar(const struct Color *restrict col1,
					      const struct Color *restrict col2,
					      int i, int j) {
	return sqrtf(SQUARE_SCALAR(col1->alpha[i] - col2->alpha[j]) +
		     SQUARE_SCALAR(col1->beta[i] - col2->beta[j]) +
		     SQUARE_SCALAR(col1->gamma[i] - col2->gamma[j]));
}

static inline float delta_cie94_diff_scalar(const struct Color *restrict colors,
					    const struct Color *restrict cents,
					    int i, int j) {

	const double K1 = 0.045;
	const double K2 = 0.015;

	const double dl = colors->alpha[i] - cents->alpha[j];
	const double da = colors->beta[i] - cents->beta[j];
	const double db = colors->gamma[i] - cents->gamma[j];
	const double C1 = sqrt(SQUARE_SCALAR(colors->beta[i]) +
			       SQUARE_SCALAR(colors->gamma[i]));
	const double C2 = sqrt(SQUARE_SCALAR(cents->beta[i]) +
			       SQUARE_SCALAR(cents->gamma[i]));
	const double dC = C1 - C2;

	const double dH2 =
	    SQUARE_SCALAR(da) + SQUARE_SCALAR(db) + SQUARE_SCALAR(dC);

	const double sl = 1.0;
	const double sc = fma(K1, C1, 1.0);
	const double sh = fma(K2, C1, 1.0);

	const double term_L = dl / sl;
	const double term_C = dC / sc;
	const double term_H = sqrt(dH2 < 0.0 ? 0.0 : dH2) / sh;

	return sqrt(SQUARE_SCALAR(term_L) + SQUARE_SCALAR(term_C) +
		    SQUARE_SCALAR(term_H));
}

static inline float
delta_ciede2000_diff_scalar(const struct Color *restrict colors,
			    const struct Color *restrict cents, int i, int j) {
	const double kl = 1;
	const double kh = 1;
	const double kc = 1;
	const double l1 = colors->alpha[i];
	const double a1 = colors->beta[i];
	const double b1 = colors->gamma[i];
	const double l2 = cents->alpha[j];
	const double a2 = cents->beta[j];
	const double b2 = cents->gamma[j];

	/* Step 1 */
	const double c1 = sqrt(SQUARE_SCALAR(a1) + SQUARE_SCALAR(b1));
	const double c2 = sqrt(SQUARE_SCALAR(a2) + SQUARE_SCALAR(b2));

	const double avgc = (c1 + c2) / 2.0;
	/* Saves calculating this twice */
	const double g =
	    0.5f * (1.0f - sqrt(POW_7_SCALAR(avgc) /
				(POW_7_SCALAR(avgc) + 6103515625.0)));
	const double a1p = (1.0 + g) * a1;
	const double a2p = (1.0 + g) * a2;
	const double c1p = sqrt(SQUARE_SCALAR(a1p) + SQUARE_SCALAR(b1));
	const double c2p = sqrt(SQUARE_SCALAR(a2p) + SQUARE_SCALAR(b2));
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
	const double dtheta =
	    30.0 * exp(-(SQUARE_SCALAR((avghp - 275.0) / 25.0)));
	const double rc =
	    2.0 * sqrt(POW_7_SCALAR(avgcp) /
		       (POW_7_SCALAR(avgcp) + 6103515625.0)); /* 6103515625.0
								  = 25^3 */
	const double sl	   = 1.0 + ((0.015 * (SQUARE_SCALAR(avglp - 50.0))) /
				    sqrt(20.0 + (SQUARE_SCALAR(avglp - 50.0))));
	const double sc	   = 1.0 + 0.045 * avgcp;
	const double sh	   = 1.0 + 0.015 * avgcp * t;
	const double rt	   = -sin(DEG2RAD(2.0 * dtheta)) * rc;
	const double terml = dlp / (kl * sl);
	const double termc = dcp / (kc * sc);
	const double termh = dHp / (kh * sh);

	return sqrt(SQUARE_SCALAR(terml) + SQUARE_SCALAR(termc) +
		    SQUARE_SCALAR(termh) + rt * termc * termh);
}

void delta_lab_euc_diff_sse(const struct Color *restrict colors,
			    const struct Color *restrict cents,
			    float *restrict diffs) {
	uint16_t num_col  = colors->num;
	uint16_t num_cent = cents->num;
	for (size_t i = 0; i + FL_PER_VEC <= num_col; i += FL_PER_VEC) {
		for (size_t j = 0; j + FL_PER_VEC <= num_cent;
		     j += FL_PER_VEC) {
			__m128 co_l = _mm_loadu_ps(&colors->alpha[i]);
			__m128 co_a = _mm_loadu_ps(&colors->beta[i]);
			__m128 co_b = _mm_loadu_ps(&colors->gamma[i]);
			__m128 ce_l = _mm_loadu_ps(&cents->alpha[j]);
			__m128 ce_a = _mm_loadu_ps(&cents->beta[j]);
			__m128 ce_b = _mm_loadu_ps(&cents->gamma[j]);

			__m128 dl	= square(_mm_sub_ps(co_l, ce_l));
			__m128 da	= square(_mm_sub_ps(co_a, ce_a));
			__m128 db	= square(_mm_sub_ps(co_b, ce_b));
			__m128 sum	= _mm_add_ps(dl, _mm_add_ps(da, db));
			__m128 diff_vec = _mm_sqrt_ps(sum);

			_mm_storeu_ps(&diffs[j * num_col + i], diff_vec);
		}
	}
	compute_tail(colors, cents, diffs, delta_lab_euc_diff_scalar);
}

void delta_cie94_diff_sse(const struct Color *restrict colors,
			  const struct Color *restrict cents,
			  float *restrict diffs) {
	const uint32_t num_col	= colors->num;
	const uint32_t num_cent = cents->num;
	const __m128   K1	= V(0.045);
	const __m128   K2	= V(0.015);
	const __m128   zero	= V(0.0f);

	for (size_t i = 0; i + FL_PER_VEC <= num_col; i += FL_PER_VEC) {
		for (size_t j = 0; j + FL_PER_VEC <= num_cent;
		     j += FL_PER_VEC) {
			__m128	     co_l = _mm_loadu_ps(&colors->alpha[i]);
			__m128	     co_a = _mm_loadu_ps(&colors->beta[i]);
			__m128	     co_b = _mm_loadu_ps(&colors->gamma[i]);
			__m128	     ce_l = _mm_loadu_ps(&cents->alpha[j]);
			__m128	     ce_a = _mm_loadu_ps(&cents->beta[j]);
			__m128	     ce_b = _mm_loadu_ps(&cents->gamma[j]);
			const __m128 dl	  = _mm_sub_ps(co_l, ce_l);
			const __m128 da	  = _mm_sub_ps(co_a, ce_a);
			const __m128 db	  = _mm_sub_ps(co_b, ce_b);
			const __m128 C1 =
			    _mm_sqrt_ps(_mm_add_ps(square(co_a), square(co_b)));
			const __m128 C2 =
			    _mm_sqrt_ps(_mm_add_ps(square(ce_a), square(ce_b)));
			const __m128 dC	 = _mm_sub_ps(C1, C2);
			const __m128 dH2 = _mm_sub_ps(
			    _mm_add_ps(square(da), square(db)), square(dC));
			const __m128 sl	    = V(1.0);
			const __m128 sc	    = FMADD(K1, C1, sl);
			const __m128 sh	    = FMADD(K2, C1, sl);
			const __m128 term_L = dl;
			const __m128 term_C = _mm_div_ps(dC, sc);
			const __m128 term_H =
			    _mm_div_ps(_mm_sqrt_ps(_mm_max_ps(dH2, zero)), sh);

			_mm_storeu_ps(
			    &diffs[j * num_col + i],
			    _mm_sqrt_ps(_mm_add_ps(
				square(term_L),
				_mm_add_ps(square(term_C), square(term_H)))));
		}
	}
	compute_tail(colors, cents, diffs, delta_cie94_diff_scalar);
}

/* The SIMD API is quite a sight, isn't it? */
void delta_ciede2000_diff_sse(const struct Color *restrict colors,
			      const struct Color *restrict cents,
			      float *restrict diffs) {
	const uint32_t num_col	= colors->num;
	const uint32_t num_cent = cents->num;
	const __m128   zero	= V(0.0f);
	const __m128   one	= V(1.0f);
	const __m128   two	= V(2.0f);
	const __m128   c180	= V(180.0f);
	const __m128   c360	= V(360.0f);

	for (size_t i = 0; i + FL_PER_VEC <= num_col; i += FL_PER_VEC) {
		for (size_t j = 0; j + FL_PER_VEC <= num_cent;
		     j += FL_PER_VEC) {
			const __m128 L1 = _mm_loadu_ps(&colors->alpha[i]);
			const __m128 a1 = _mm_loadu_ps(&colors->beta[i]);
			const __m128 b1 = _mm_loadu_ps(&colors->gamma[i]);
			const __m128 L2 = _mm_loadu_ps(&cents->alpha[j]);
			const __m128 a2 = _mm_loadu_ps(&cents->beta[j]);
			const __m128 b2 = _mm_loadu_ps(&cents->gamma[j]);

			/* Step 1: Find Cip, hip*/
			const __m128 C1 =
			    _mm_sqrt_ps(_mm_add_ps(square(a1), square(b1)));
			const __m128 C2 =
			    _mm_sqrt_ps(_mm_add_ps(square(a2), square(b2)));
			const __m128 avgC = _mm_div_ps(_mm_add_ps(C1, C2), two);
			/* This is gross */
			const __m128 g = _mm_mul_ps(
			    V(0.5f),
			    _mm_sub_ps(one,
				       _mm_sqrt_ps(_mm_div_ps(
					   pow7(avgC),
					   _mm_add_ps(pow7(avgC),
						      V(6103515625.0f))))));
			const __m128 a1p = _mm_mul_ps(_mm_add_ps(one, g), a1);
			const __m128 a2p = _mm_mul_ps(_mm_add_ps(one, g), a2);
			const __m128 C1p =
			    _mm_sqrt_ps(_mm_add_ps(square(a1p), square(b1)));
			const __m128 C2p =
			    _mm_sqrt_ps(_mm_add_ps(square(a2p), square(b2)));
			__m128	     mask = _mm_and_ps(_mm_cmpeq_ps(b1, zero),
						       _mm_cmpeq_ps(a1p, zero));
			const __m128 h1p  = _mm_blendv_ps(
			     rad_2_deg(_mm_atan2_ps(b1, a1p)), zero, mask);
			mask		 = _mm_and_ps(_mm_cmpeq_ps(a2p, zero),
						      _mm_cmpeq_ps(b2, zero));
			const __m128 h2p = _mm_blendv_ps(
			    rad_2_deg(_mm_atan2_ps(b2, a2p)), zero, mask);

			/* Step 2: Find dLp, dCp, dHp */
			const __m128 dLp = _mm_sub_ps(L2, L1);
			const __m128 dCp = _mm_sub_ps(C2p, C1p);
			__m128	     dhp = _mm_setzero_ps();
			mask		 = _mm_and_ps(
				_mm_cmpneq_ps(_mm_mul_ps(C1p, C2p), zero),
				_mm_cmple_ps(_mm_abs_ps(_mm_sub_ps(h2p, h1p)),
						     c180));
			dhp  = _mm_blendv_ps(dhp, _mm_sub_ps(h2p, h1p), mask);
			mask = _mm_and_ps(
			    _mm_cmpneq_ps(_mm_mul_ps(C1p, C2p), zero),
			    _mm_cmpgt_ps(_mm_sub_ps(h2p, h1p), c180));
			dhp = _mm_blendv_ps(
			    dhp, _mm_sub_ps(h2p, _mm_sub_ps(h1p, c360)), mask);
			mask = _mm_and_ps(
			    _mm_cmpneq_ps(_mm_mul_ps(C1p, C2p), zero),
			    _mm_cmplt_ps(_mm_sub_ps(h2p, h1p),
					 _mm_mul_ps(V(-1.0f), c180)));
			dhp = _mm_blendv_ps(
			    dhp, _mm_sub_ps(h2p, _mm_add_ps(h1p, c360)), mask);
			const __m128 dHp = _mm_mul_ps(
			    _mm_mul_ps(two, _mm_sqrt_ps(_mm_mul_ps(C1p, C2p))),
			    _mm_sin_ps(_mm_div_ps(dhp, two)));

			/* Step 3: Fin */
			const __m128 avgLp =
			    _mm_div_ps(_mm_add_ps(L1, L2), two);
			const __m128 avgCp =
			    _mm_div_ps(_mm_add_ps(C1p, C2p), two);

			__m128 avghp = _mm_add_ps(h1p, h2p);
			mask	     = _mm_and_ps(
			    _mm_cmple_ps(_mm_abs_ps(_mm_sub_ps(h1p, h2p)),
						 c180),
			    _mm_cmpneq_ps(_mm_mul_ps(C1p, C2p), zero));
			avghp = _mm_blendv_ps(
			    avghp, _mm_div_ps(_mm_add_ps(h1p, h2p), two), mask);
			mask = _mm_and_ps(
			    _mm_cmpgt_ps(_mm_abs_ps(_mm_sub_ps(h1p, h2p)),
					 c180),
			    _mm_cmplt_ps(_mm_add_ps(h1p, h2p), c360));
			avghp = _mm_blendv_ps(
			    avghp,
			    _mm_div_ps(_mm_add_ps(h1p, _mm_add_ps(h2p, c360)),
				       two),
			    mask);
			mask = _mm_and_ps(
			    _mm_cmpge_ps(_mm_abs_ps(_mm_sub_ps(h1p, h2p)),
					 c180),
			    _mm_and_ps(
				_mm_cmpge_ps(_mm_add_ps(h1p, h2p), c360),
				_mm_cmpeq_ps(_mm_mul_ps(C1p, C2p), zero)));
			avghp = _mm_blendv_ps(
			    avghp,
			    _mm_div_ps(_mm_sub_ps(_mm_add_ps(h1p, h2p), c360),
				       two),
			    mask);

			__m128 T = _mm_sub_ps(
			    one, _mm_mul_ps(V(0.17f),
					    _mm_cos_ps(deg_2_rad(
						_mm_sub_ps(avghp, V(30.0f))))));
			T = _mm_add_ps(
			    T,
			    _mm_mul_ps(V(0.24f), _mm_cos_ps(deg_2_rad(
						     _mm_mul_ps(two, avghp)))));
			T = _mm_add_ps(
			    T, _mm_mul_ps(
				   V(0.32f),
				   _mm_cos_ps(deg_2_rad(_mm_add_ps(
				       _mm_mul_ps(V(3.0f), avghp), V(6.0f))))));
			T = _mm_sub_ps(
			    T, _mm_mul_ps(V(20.0f),
					  _mm_cos_ps(deg_2_rad(_mm_sub_ps(
					      _mm_mul_ps(V(4.0f), avghp),
					      V(63.0f))))));
			const __m128 dtheta = _mm_mul_ps(
			    V(30.0f),
			    _mm_exp_ps(_mm_mul_ps(
				V(-1.0f),
				square(_mm_div_ps(_mm_sub_ps(avghp, V(275.0f)),
						  V(25.0f))))));
			const __m128 RC = _mm_mul_ps(
			    V(2.0f),
			    _mm_sqrt_ps(_mm_div_ps(
				pow7(avgCp),
				_mm_add_ps(pow7(avgCp), V(6103515625.0f)))));
			const __m128 SL = _mm_add_ps(
			    one,
			    _mm_div_ps(
				_mm_mul_ps(V(0.015f),
					   square(_mm_sub_ps(avgLp, V(50.0f)))),
				_mm_sqrt_ps(_mm_add_ps(
				    V(20.0f),
				    square(_mm_sub_ps(avgLp, V(50.0f)))))));
			const __m128 SC =
			    _mm_add_ps(one, _mm_mul_ps(V(0.045f), avgCp));
			const __m128 SH = _mm_add_ps(
			    one, _mm_mul_ps(V(0.015), _mm_mul_ps(avgCp, T)));
			const __m128 RT = _mm_mul_ps(
			    V(-1.0f), _mm_mul_ps(_mm_sin_ps(deg_2_rad(
						     _mm_mul_ps(two, dtheta))),
						 RC));
			/* Omit KL, KH, KC as they are all just 1 for our
			 * purposes*/
			__m128 diff = square(_mm_div_ps(dLp, SL));
			diff = _mm_add_ps(diff, square(_mm_div_ps(dCp, SC)));
			diff = _mm_add_ps(diff, square(_mm_div_ps(dHp, SH)));
			diff = _mm_add_ps(
			    diff,
			    _mm_mul_ps(RT, _mm_mul_ps(_mm_div_ps(dCp, SC),
						      _mm_div_ps(dHp, SH))));
			diff = _mm_sqrt_ps(diff);
			_mm_storeu_ps(&diffs[j * num_col + i], diff);
		}
	}
	compute_tail(colors, cents, diffs, delta_ciede2000_diff_scalar);
}
