#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <immintrin.h>
#include <math.h>

#include "Color.h"

#define SQUARE_SCALAR(x) ((x) * (x))
#define POW_7_SCALAR(x) ((x) * (x) * (x) * (x) * (x) * (x) * (x))
#define DEG2RAD(x) (((x) / (180.0)) * (M_PI))
#define RAD2DEG(x) (((x) / M_PI) * (180.0))
#define _mm256_abs_ps(x) _mm256_and_ps((x), _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF)))
#define FL_PER_V256EC 8
#define _MM256_PI 3.14159265f
#define _MM256_2PI 6.28318531f
#define _MM256_PI_2 1.57079633f
#define _MM256_1_OV256ER_2PI 0.159154943f

static inline __m256 pow7(__m256 x) {
	__m256 x2 = _mm256_mul_ps(x, x);
	__m256 x4 = _mm256_mul_ps(x2, x2);
	return _mm256_mul_ps(x, _mm256_mul_ps(x2, x4));
}

static inline __m256 square(__m256 x) {
	return _mm256_mul_ps(x, x);
}

static inline __m256 deg_2_rad(__m256 x) {
	return _mm256_mul_ps(_mm256_div_ps(x, _mm256_set1_ps(180.0f)),
			     _mm256_set1_ps(_MM256_PI));
}

static inline __m256 rad_2_deg(__m256 x) {
	return _mm256_mul_ps(_mm256_div_ps(x, _mm256_set1_ps(_MM256_PI)),
			     _mm256_set1_ps(180.0f));
}

static inline __m256 V256(float x) {
	return _mm256_set1_ps(x);
}

static inline void
compute_tail(const struct Color *restrict colors,
	     const struct Color *restrict cents, float *diffs,
	     float (*scalar_diff)(const struct Color *restrict,
				  const struct Color *restrict, int, int)) {
	const uint32_t num_col	     = colors->num;
	const uint32_t num_cent	     = cents->num;
	size_t	       col_not_comp  = num_col - (num_col % FL_PER_V256EC);
	size_t	       cent_not_comp = num_cent - (num_cent % FL_PER_V256EC);
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

static inline __m256 _mm256_atan_ps(__m256 x) {
	const __m256 one  = V256(1.0f);
	const __m256 pi_2 = V256(_MM256_PI_2);

	__m256 sign_mask = _mm256_and_ps(x, V256(-0.0f));
	__m256 abs_x	 = _mm256_andnot_ps(sign_mask, x);

	__m256 use_identity_mask = _mm256_cmp_ps(abs_x, one, _CMP_GT_OQ);

	__m256 x_inv = _mm256_div_ps(one, abs_x);

	__m256 x_red = _mm256_or_ps(_mm256_and_ps(use_identity_mask, x_inv),
				    _mm256_andnot_ps(use_identity_mask, abs_x));

	__m256 x2 = _mm256_mul_ps(x_red, x_red);
	__m256 y  = V256(-0.01891985f);
	y	  = _mm256_fmadd_ps(y, x2, V256(0.06909289f));
	y	  = _mm256_fmadd_ps(y, x2, V256(-0.12945813f));
	y	  = _mm256_fmadd_ps(y, x2, V256(0.19787976f));
	y	  = _mm256_fmadd_ps(y, x2, V256(-0.33319414f));
	y	  = _mm256_fmadd_ps(y, x2, V256(0.99999750f));

	y = _mm256_mul_ps(y, x_red);

	__m256 atan_corrected = _mm256_sub_ps(pi_2, y);
	y = _mm256_or_ps(_mm256_and_ps(use_identity_mask, atan_corrected),
			 _mm256_andnot_ps(use_identity_mask, y));

	y = _mm256_or_ps(y, sign_mask);
	return y;
}

static inline __m256 _mm256_atan2_ps(__m256 y, __m256 x) {
	const __m256 zero     = V256(0.0f);
	const __m256 pi	      = V256(_MM256_PI);
	const __m256 pi_2     = V256(_MM256_PI_2);
	const __m256 neg_pi_2 = V256(-_MM256_PI_2);

	__m256 abs_x = _mm256_andnot_ps(V256(-0.0f), x);

	__m256 safe_x = _mm256_blendv_ps(
	    x, V256(1.0f), _mm256_cmp_ps(abs_x, zero, _CMP_EQ_OQ));
	__m256 z      = _mm256_div_ps(y, safe_x);
	__m256 atan_z = _mm256_atan_ps(z);

	__m256 x_lt0 = _mm256_cmp_ps(x, zero, _CMP_LT_OQ);
	__m256 y_eq0 = _mm256_cmp_ps(y, zero, _CMP_EQ_OQ);
	__m256 y_lt0 = _mm256_cmp_ps(y, zero, _CMP_LT_OQ);
	__m256 y_gt0 = _mm256_cmp_ps(y, zero, _CMP_GT_OQ);
	__m256 x_eq0 = _mm256_cmp_ps(x, zero, _CMP_EQ_OQ);

	__m256 add_pi = _mm256_and_ps(x_lt0, y_gt0);
	__m256 sub_pi = _mm256_and_ps(x_lt0, y_lt0);
	__m256 result = atan_z;
	result	      = _mm256_add_ps(result, _mm256_and_ps(add_pi, pi));
	result	      = _mm256_sub_ps(result, _mm256_and_ps(sub_pi, pi));

	__m256 x0_y0   = _mm256_and_ps(x_eq0, y_eq0);
	__m256 x0_ygt0 = _mm256_and_ps(x_eq0, y_gt0);
	__m256 x0_ylt0 = _mm256_and_ps(x_eq0, y_lt0);

	result = _mm256_blendv_ps(result, zero, x0_y0);
	result = _mm256_blendv_ps(result, pi_2, x0_ygt0);
	result = _mm256_blendv_ps(result, neg_pi_2, x0_ylt0);

	return result;
}

static inline __m256 _mm256_sin_ps(__m256 x) {
	const __m256 inv_2pi = V256(_MM256_1_OV256ER_2PI);
	const __m256 two_pi  = V256(_MM256_2PI);
	const __m256 pi	     = V256(_MM256_PI);

	__m256 k = _mm256_mul_ps(x, inv_2pi);
	k	 = _mm256_round_ps(k, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
	x	 = _mm256_sub_ps(x, _mm256_mul_ps(k, two_pi));
	x	 = _mm256_sub_ps(x, pi);

	__m256 x2 = _mm256_mul_ps(x, x);
	__m256 y  = V256(0.00000227f);
	y	  = _mm256_fmadd_ps(y, x2, V256(-0.00019506f));
	y	  = _mm256_fmadd_ps(y, x2, V256(0.00832404f));
	y	  = _mm256_fmadd_ps(y, x2, V256(-0.16665670f));
	y	  = _mm256_fmadd_ps(y, x2, V256(0.99999708f));

	y = _mm256_mul_ps(y, -x);
	return _mm256_min_ps(_mm256_max_ps(y, V256(-1.0f)), V256(1.0f));
}

static inline __m256 _mm256_cos_ps(__m256 x) {
	const __m256 inv_2pi = V256(_MM256_1_OV256ER_2PI);
	const __m256 two_pi  = V256(_MM256_2PI);
	const __m256 pi	     = V256(_MM256_PI);

	__m256 k = _mm256_mul_ps(x, inv_2pi);
	k	 = _mm256_round_ps(k, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
	x	 = _mm256_sub_ps(x, _mm256_mul_ps(k, two_pi));
	x	 = _mm256_sub_ps(x, pi);

	__m256 x2 = _mm256_mul_ps(x, x);
	__m256 y  = V256(0.00001974f);
	y	  = _mm256_fmadd_ps(y, x2, V256(-0.00135731f));
	y	  = _mm256_fmadd_ps(y, x2, V256(0.04159290f));
	y	  = _mm256_fmadd_ps(y, x2, V256(-0.49994522f));
	y	  = _mm256_fmadd_ps(y, x2, V256(0.99999571f));

	y = _mm256_mul_ps(y, V256(-1.0f));
	return _mm256_min_ps(_mm256_max_ps(y, V256(-1.0f)), V256(1.0f));
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
static inline __m256 _mm256_exp_ps(__m256 x) {
	x	  = _mm256_min_ps(_mm256_max_ps(x, V256(-87.33654)), V256(88.72283));
	__m256	a = V256(12102203.0f); /* (1 << 23) / log(2) */
	__m256i b = _mm256_set1_epi32(127 * (1 << 23) - 298765);
	__m256i t =
	    _mm256_add_epi32(_mm256_cvtps_epi32(_mm256_mul_ps(a, x)), b);
	return _mm256_castsi256_ps(t);
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

void delta_lab_euc_diff_avx(const struct Color *restrict colors,
			    const struct Color *restrict cents,
			    float *restrict diffs) {
	uint16_t num_col  = colors->num;
	uint16_t num_cent = cents->num;
	for (size_t i = 0; i + FL_PER_V256EC <= num_col; i += FL_PER_V256EC) {
		for (size_t j = 0; j + FL_PER_V256EC <= num_cent;
		     j += FL_PER_V256EC) {
			__m256 co_l = _mm256_loadu_ps(&colors->alpha[i]);
			__m256 co_a = _mm256_loadu_ps(&colors->beta[i]);
			__m256 co_b = _mm256_loadu_ps(&colors->gamma[i]);
			__m256 ce_l = _mm256_loadu_ps(&cents->alpha[j]);
			__m256 ce_a = _mm256_loadu_ps(&cents->beta[j]);
			__m256 ce_b = _mm256_loadu_ps(&cents->gamma[j]);

			__m256 dl  = square(_mm256_sub_ps(co_l, ce_l));
			__m256 da  = square(_mm256_sub_ps(co_a, ce_a));
			__m256 db  = square(_mm256_sub_ps(co_b, ce_b));
			__m256 sum = _mm256_add_ps(dl, _mm256_add_ps(da, db));
			__m256 diff_vec = _mm256_sqrt_ps(sum);

			_mm256_storeu_ps(&diffs[j * num_col + i], diff_vec);
		}
	}
	compute_tail(colors, cents, diffs, delta_lab_euc_diff_scalar);
}

void delta_cie94_diff_avx(const struct Color *restrict colors,
			  const struct Color *restrict cents,
			  float *restrict diffs) {
	const uint32_t num_col	= colors->num;
	const uint32_t num_cent = cents->num;
	const __m256   K1	= V256(0.045);
	const __m256   K2	= V256(0.015);
	const __m256   zero	= V256(0.0f);

	for (size_t i = 0; i + FL_PER_V256EC <= num_col; i += FL_PER_V256EC) {
		for (size_t j = 0; j + FL_PER_V256EC <= num_cent;
		     j += FL_PER_V256EC) {
			__m256	     co_l = _mm256_loadu_ps(&colors->alpha[i]);
			__m256	     co_a = _mm256_loadu_ps(&colors->beta[i]);
			__m256	     co_b = _mm256_loadu_ps(&colors->gamma[i]);
			__m256	     ce_l = _mm256_loadu_ps(&cents->alpha[j]);
			__m256	     ce_a = _mm256_loadu_ps(&cents->beta[j]);
			__m256	     ce_b = _mm256_loadu_ps(&cents->gamma[j]);
			
			const __m256 dl	  = _mm256_sub_ps(co_l, ce_l);
			const __m256 da	  = _mm256_sub_ps(co_a, ce_a);
			const __m256 db	  = _mm256_sub_ps(co_b, ce_b);
			const __m256 C1	  = _mm256_sqrt_ps(
			      _mm256_add_ps(square(co_a), square(co_b)));
			const __m256 C2 = _mm256_sqrt_ps(
			    _mm256_add_ps(square(ce_a), square(ce_b)));
			const __m256 dC	 = _mm256_sub_ps(C1, C2);
			const __m256 dH2 = _mm256_sub_ps(
			    _mm256_add_ps(square(da), square(db)), square(dC));
			const __m256 sl	    = V256(1.0);
			const __m256 sc	    = _mm256_fmadd_ps(K1, C1, sl);
			const __m256 sh	    = _mm256_fmadd_ps(K2, C1, sl);
			const __m256 term_L = dl;
			const __m256 term_C = _mm256_div_ps(dC, sc);
			const __m256 term_H = _mm256_div_ps(
			    _mm256_sqrt_ps(_mm256_max_ps(dH2, zero)), sh);

			_mm256_storeu_ps(&diffs[j * num_col + i],
					 _mm256_sqrt_ps(_mm256_add_ps(
					     square(term_L),
					     _mm256_add_ps(square(term_C),
							   square(term_H)))));
		}
	}
	compute_tail(colors, cents, diffs, delta_cie94_diff_scalar);
}

/* The SIMD API is quite a sight, isn't it? */
void delta_ciede2000_diff_avx(const struct Color *restrict colors,
			      const struct Color *restrict cents,
			      float *restrict diffs) {
	const uint32_t num_col	= colors->num;
	const uint32_t num_cent = cents->num;
	const __m256   zero	= V256(0.0f);
	const __m256   one	= V256(1.0f);
	const __m256   two	= V256(2.0f);
	const __m256   c180	= V256(180.0f);
	const __m256   c360	= V256(360.0f);

	for (size_t i = 0; i + FL_PER_V256EC <= num_col; i += FL_PER_V256EC) {
		for (size_t j = 0; j + FL_PER_V256EC <= num_cent;
		     j += FL_PER_V256EC) {
			const __m256 L1 = _mm256_loadu_ps(&colors->alpha[i]);
			const __m256 a1 = _mm256_loadu_ps(&colors->beta[i]);
			const __m256 b1 = _mm256_loadu_ps(&colors->gamma[i]);
			const __m256 L2 = _mm256_loadu_ps(&cents->alpha[j]);
			const __m256 a2 = _mm256_loadu_ps(&cents->beta[j]);
			const __m256 b2 = _mm256_loadu_ps(&cents->gamma[j]);

			/* Step 1: Find Cip, hip*/
			const __m256 C1 = _mm256_sqrt_ps(
			    _mm256_add_ps(square(a1), square(b1)));
			const __m256 C2 = _mm256_sqrt_ps(
			    _mm256_add_ps(square(a2), square(b2)));
			const __m256 avgC =
			    _mm256_div_ps(_mm256_add_ps(C1, C2), two);
			/* This is gross */
			const __m256 g = _mm256_mul_ps(
			    V256(0.5f),
			    _mm256_sub_ps(
				one, _mm256_sqrt_ps(_mm256_div_ps(
					 pow7(avgC),
					 _mm256_add_ps(pow7(avgC),
						       V256(6103515625.0f))))));
			const __m256 a1p =
			    _mm256_mul_ps(_mm256_add_ps(one, g), a1);
			const __m256 a2p =
			    _mm256_mul_ps(_mm256_add_ps(one, g), a2);
			const __m256 C1p = _mm256_sqrt_ps(
			    _mm256_add_ps(square(a1p), square(b1)));
			const __m256 C2p = _mm256_sqrt_ps(
			    _mm256_add_ps(square(a2p), square(b2)));
			__m256 mask =
			    _mm256_and_ps(_mm256_cmp_ps(b1, zero, _CMP_EQ_OQ),
					  _mm256_cmp_ps(a1p, zero, _CMP_EQ_OQ));
			const __m256 h1p = _mm256_blendv_ps(
			    rad_2_deg(_mm256_atan2_ps(b1, a1p)), zero, mask);
			mask =
			    _mm256_and_ps(_mm256_cmp_ps(a2p, zero, _CMP_EQ_OQ),
					  _mm256_cmp_ps(b2, zero, _CMP_EQ_OQ));
			const __m256 h2p = _mm256_blendv_ps(
			    rad_2_deg(_mm256_atan2_ps(b2, a2p)), zero, mask);

			/* Step 2: Find dLp, dCp, dHp */
			const __m256 dLp = _mm256_sub_ps(L2, L1);
			const __m256 dCp = _mm256_sub_ps(C2p, C1p);
			__m256	     dhp = _mm256_setzero_ps();
			mask		 = _mm256_and_ps(
				_mm256_cmp_ps(_mm256_mul_ps(C1p, C2p), zero,
						      _CMP_NEQ_OQ),
				_mm256_cmp_ps(
				    _mm256_abs_ps(_mm256_sub_ps(h2p, h1p)), c180,
				    _CMP_LE_OQ));
			dhp = _mm256_blendv_ps(dhp, _mm256_sub_ps(h2p, h1p),
					       mask);
			mask =
			    _mm256_and_ps(_mm256_cmp_ps(_mm256_mul_ps(C1p, C2p),
							zero, _CMP_NEQ_OQ),
					  _mm256_cmp_ps(_mm256_sub_ps(h2p, h1p),
							c180, _CMP_GT_OQ));
			dhp = _mm256_blendv_ps(
			    dhp, _mm256_sub_ps(h2p, _mm256_sub_ps(h1p, c360)),
			    mask);
			mask = _mm256_and_ps(
			    _mm256_cmp_ps(_mm256_mul_ps(C1p, C2p), zero,
					  _CMP_NEQ_OQ),
			    _mm256_cmp_ps(_mm256_sub_ps(h2p, h1p),
					  _mm256_mul_ps(V256(-1.0f), c180),
					  _CMP_LT_OQ));
			dhp = _mm256_blendv_ps(
			    dhp, _mm256_sub_ps(h2p, _mm256_add_ps(h1p, c360)),
			    mask);
			const __m256 dHp = _mm256_mul_ps(
			    _mm256_mul_ps(
				two, _mm256_sqrt_ps(_mm256_mul_ps(C1p, C2p))),
			    _mm256_sin_ps(_mm256_div_ps(dhp, two)));

			/* Step 3: Fin */
			const __m256 avgLp =
			    _mm256_div_ps(_mm256_add_ps(L1, L2), two);
			const __m256 avgCp =
			    _mm256_div_ps(_mm256_add_ps(C1p, C2p), two);

			__m256 avghp = _mm256_add_ps(h1p, h2p);
			mask	     = _mm256_and_ps(
			    _mm256_cmp_ps(
				_mm256_abs_ps(_mm256_sub_ps(h1p, h2p)), c180,
				_CMP_LE_OQ),
			    _mm256_cmp_ps(_mm256_mul_ps(C1p, C2p), zero,
						  _CMP_NEQ_OQ));
			avghp = _mm256_blendv_ps(
			    avghp, _mm256_div_ps(_mm256_add_ps(h1p, h2p), two),
			    mask);
			mask = _mm256_and_ps(
			    _mm256_cmp_ps(
				_mm256_abs_ps(_mm256_sub_ps(h1p, h2p)), c180,
				_CMP_GT_OQ),
			    _mm256_cmp_ps(_mm256_add_ps(h1p, h2p), c360,
					  _CMP_LT_OQ));
			avghp = _mm256_blendv_ps(
			    avghp,
			    _mm256_div_ps(
				_mm256_add_ps(h1p, _mm256_add_ps(h2p, c360)),
				two),
			    mask);
			mask = _mm256_and_ps(
			    _mm256_cmp_ps(
				_mm256_abs_ps(_mm256_sub_ps(h1p, h2p)), c180,
				_CMP_GE_OQ),
			    _mm256_and_ps(_mm256_cmp_ps(_mm256_add_ps(h1p, h2p),
							c360, _CMP_GE_OQ),
					  _mm256_cmp_ps(_mm256_mul_ps(C1p, C2p),
							zero, _CMP_EQ_OQ)));
			avghp = _mm256_blendv_ps(
			    avghp,
			    _mm256_div_ps(
				_mm256_sub_ps(_mm256_add_ps(h1p, h2p), c360),
				two),
			    mask);

			__m256 T = _mm256_sub_ps(
			    one,
			    _mm256_mul_ps(V256(0.17f),
					  _mm256_cos_ps(deg_2_rad(_mm256_sub_ps(
					      avghp, V256(30.0f))))));
			T = _mm256_add_ps(
			    T, _mm256_mul_ps(V256(0.24f),
					     _mm256_cos_ps(deg_2_rad(
						 _mm256_mul_ps(two, avghp)))));
			T = _mm256_add_ps(
			    T,
			    _mm256_mul_ps(
				V256(0.32f),
				_mm256_cos_ps(deg_2_rad(_mm256_add_ps(
				    _mm256_mul_ps(V256(3.0f), avghp), V256(6.0f))))));
			T = _mm256_sub_ps(
			    T,
			    _mm256_mul_ps(V256(20.0f),
					  _mm256_cos_ps(deg_2_rad(_mm256_sub_ps(
					      _mm256_mul_ps(V256(4.0f), avghp),
					      V256(63.0f))))));
			const __m256 dtheta = _mm256_mul_ps(
			    V256(30.0f),
			    _mm256_exp_ps(_mm256_mul_ps(
				V256(-1.0f), square(_mm256_div_ps(
					      _mm256_sub_ps(avghp, V256(275.0f)),
					      V256(25.0f))))));
			const __m256 RC = _mm256_mul_ps(
			    V256(2.0f),
			    _mm256_sqrt_ps(_mm256_div_ps(
				pow7(avgCp),
				_mm256_add_ps(pow7(avgCp), V256(6103515625.0f)))));
			const __m256 SL = _mm256_add_ps(
			    one,
			    _mm256_div_ps(
				_mm256_mul_ps(V256(0.015f), square(_mm256_sub_ps(
							     avgLp, V256(50.0f)))),
				_mm256_sqrt_ps(_mm256_add_ps(
				    V256(20.0f),
				    square(_mm256_sub_ps(avgLp, V256(50.0f)))))));
			const __m256 SC =
			    _mm256_add_ps(one, _mm256_mul_ps(V256(0.045f), avgCp));
			const __m256 SH = _mm256_add_ps(
			    one,
			    _mm256_mul_ps(V256(0.015), _mm256_mul_ps(avgCp, T)));
			const __m256 RT = _mm256_mul_ps(
			    V256(-1.0f),
			    _mm256_mul_ps(_mm256_sin_ps(deg_2_rad(
					      _mm256_mul_ps(two, dtheta))),
					  RC));
			/* Omit KL, KH, KC as they are all just 1 for our
			 * purposes*/
			__m256 diff = square(_mm256_div_ps(dLp, SL));
			diff =
			    _mm256_add_ps(diff, square(_mm256_div_ps(dCp, SC)));
			diff =
			    _mm256_add_ps(diff, square(_mm256_div_ps(dHp, SH)));
			diff = _mm256_add_ps(
			    diff,
			    _mm256_mul_ps(
				RT, _mm256_mul_ps(_mm256_div_ps(dCp, SC),
						  _mm256_div_ps(dHp, SH))));
			diff = _mm256_sqrt_ps(diff);
			_mm256_storeu_ps(&diffs[j * num_col + i], diff);
		}
	}
	compute_tail(colors, cents, diffs, delta_ciede2000_diff_scalar);
}
