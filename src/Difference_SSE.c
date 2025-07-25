#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "Color.h"
#include "Math_SSE.h"

#define SQUARE_SCALAR(x) ((x) * (x))
#define POW_7_SCALAR(x) ((x) * (x) * (x) * (x) * (x) * (x) * (x))
#define DEG2RAD(x) (((x) / (180.0)) * (M_PI))
#define RAD2DEG(x) (((x) / M_PI) * (180.0))
#define FL_PER_VEC 4

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

			__m128 dl	= _mm_square_ps(_mm_sub_ps(co_l, ce_l));
			__m128 da	= _mm_square_ps(_mm_sub_ps(co_a, ce_a));
			__m128 db	= _mm_square_ps(_mm_sub_ps(co_b, ce_b));
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
			const __m128 C1	  = _mm_sqrt_ps(_mm_add_ps(
			      _mm_square_ps(co_a), _mm_square_ps(co_b)));
			const __m128 C2	  = _mm_sqrt_ps(_mm_add_ps(
			      _mm_square_ps(ce_a), _mm_square_ps(ce_b)));
			const __m128 dC	  = _mm_sub_ps(C1, C2);
			const __m128 dH2  = _mm_sub_ps(
			     _mm_add_ps(_mm_square_ps(da), _mm_square_ps(db)),
			     _mm_square_ps(dC));
			const __m128 sl	    = V(1.0);
			const __m128 sc	    = FMADD(K1, C1, sl);
			const __m128 sh	    = FMADD(K2, C1, sl);
			const __m128 term_L = dl;
			const __m128 term_C = _mm_div_ps(dC, sc);
			const __m128 term_H =
			    _mm_div_ps(_mm_sqrt_ps(_mm_max_ps(dH2, zero)), sh);

			_mm_storeu_ps(&diffs[j * num_col + i],
				      _mm_sqrt_ps(_mm_add_ps(
					  _mm_square_ps(term_L),
					  _mm_add_ps(_mm_square_ps(term_C),
						     _mm_square_ps(term_H)))));
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
			const __m128 C1 = _mm_sqrt_ps(
			    _mm_add_ps(_mm_square_ps(a1), _mm_square_ps(b1)));
			const __m128 C2 = _mm_sqrt_ps(
			    _mm_add_ps(_mm_square_ps(a2), _mm_square_ps(b2)));
			const __m128 avgC = _mm_div_ps(_mm_add_ps(C1, C2), two);
			/* This is gross */
			const __m128 g = _mm_mul_ps(
			    V(0.5f),
			    _mm_sub_ps(one,
				       _mm_sqrt_ps(_mm_div_ps(
					   _mm_pow7_ps(avgC),
					   _mm_add_ps(_mm_pow7_ps(avgC),
						      V(6103515625.0f))))));
			const __m128 a1p = _mm_mul_ps(_mm_add_ps(one, g), a1);
			const __m128 a2p = _mm_mul_ps(_mm_add_ps(one, g), a2);
			const __m128 C1p = _mm_sqrt_ps(
			    _mm_add_ps(_mm_square_ps(a1p), _mm_square_ps(b1)));
			const __m128 C2p = _mm_sqrt_ps(
			    _mm_add_ps(_mm_square_ps(a2p), _mm_square_ps(b2)));
			__m128	     mask = _mm_and_ps(_mm_cmpeq_ps(b1, zero),
						       _mm_cmpeq_ps(a1p, zero));
			const __m128 h1p  = _mm_blendv_ps(
			     _mm_rad_2_deg_ps(_mm_atan2_ps(b1, a1p)), zero,
			     mask);
			mask		 = _mm_and_ps(_mm_cmpeq_ps(a2p, zero),
						      _mm_cmpeq_ps(b2, zero));
			const __m128 h2p = _mm_blendv_ps(
			    _mm_rad_2_deg_ps(_mm_atan2_ps(b2, a2p)), zero,
			    mask);

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
					    _mm_cos_ps(_mm_deg_2_rad_ps(
						_mm_sub_ps(avghp, V(30.0f))))));
			T = _mm_add_ps(
			    T,
			    _mm_mul_ps(V(0.24f), _mm_cos_ps(_mm_deg_2_rad_ps(
						     _mm_mul_ps(two, avghp)))));
			T = _mm_add_ps(
			    T, _mm_mul_ps(
				   V(0.32f),
				   _mm_cos_ps(_mm_deg_2_rad_ps(_mm_add_ps(
				       _mm_mul_ps(V(3.0f), avghp), V(6.0f))))));
			T = _mm_sub_ps(
			    T,
			    _mm_mul_ps(
				V(20.0f),
				_mm_cos_ps(_mm_deg_2_rad_ps(_mm_sub_ps(
				    _mm_mul_ps(V(4.0f), avghp), V(63.0f))))));
			const __m128 dtheta = _mm_mul_ps(
			    V(30.0f),
			    _mm_exp_ps(_mm_mul_ps(
				V(-1.0f),
				_mm_square_ps(_mm_div_ps(
				    _mm_sub_ps(avghp, V(275.0f)), V(25.0f))))));
			const __m128 RC = _mm_mul_ps(
			    V(2.0f), _mm_sqrt_ps(_mm_div_ps(
					 _mm_pow7_ps(avgCp),
					 _mm_add_ps(_mm_pow7_ps(avgCp),
						    V(6103515625.0f)))));
			const __m128 SL = _mm_add_ps(
			    one,
			    _mm_div_ps(
				_mm_mul_ps(V(0.015f), _mm_square_ps(_mm_sub_ps(
							  avgLp, V(50.0f)))),
				_mm_sqrt_ps(_mm_add_ps(
				    V(20.0f), _mm_square_ps(_mm_sub_ps(
						  avgLp, V(50.0f)))))));
			const __m128 SC =
			    _mm_add_ps(one, _mm_mul_ps(V(0.045f), avgCp));
			const __m128 SH = _mm_add_ps(
			    one, _mm_mul_ps(V(0.015), _mm_mul_ps(avgCp, T)));
			const __m128 RT = _mm_mul_ps(
			    V(-1.0f), _mm_mul_ps(_mm_sin_ps(_mm_deg_2_rad_ps(
						     _mm_mul_ps(two, dtheta))),
						 RC));
			/* Omit KL, KH, KC as they are all just 1 for our
			 * purposes*/
			__m128 diff = _mm_square_ps(_mm_div_ps(dLp, SL));
			diff	    = _mm_add_ps(diff,
						 _mm_square_ps(_mm_div_ps(dCp, SC)));
			diff	    = _mm_add_ps(diff,
						 _mm_square_ps(_mm_div_ps(dHp, SH)));
			diff	    = _mm_add_ps(
				   diff,
				   _mm_mul_ps(RT, _mm_mul_ps(_mm_div_ps(dCp, SC),
							     _mm_div_ps(dHp, SH))));
			diff = _mm_sqrt_ps(diff);
			_mm_storeu_ps(&diffs[j * num_col + i], diff);
		}
	}
	compute_tail(colors, cents, diffs, delta_ciede2000_diff_scalar);
}
