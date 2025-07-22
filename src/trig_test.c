#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#if defined(__x86_64__)
#include <x86intrin.h>
#endif
#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

#define DELTA 1e-1
#define _MM_PI 3.14159265f
#define _MM_2PI 6.28318531f
#define _MM_PI_2 1.57079633f
#define _MM_1_OVER_2PI 0.159154943f

static const uint64_t RUNS = 134217728;

uint64_t time_trig(float (*func)(float), const float *x,
		   const uint64_t num_runs) {
	/* Warm the function up */
	for (int i = 0; i < 1000; ++i) {
		func(x[i]);
	}

	uint64_t tic = __rdtsc();
	for (uint64_t i = 0; i < num_runs; ++i) {

		func(x[i]);
	}
	uint64_t toc = __rdtsc();

	return (toc - tic);
}

uint64_t time_trig_sse(__m128 (*func)(__m128), const float *x,
		       const uint64_t num_runs) {
	uint64_t tic = __rdtsc();

	for (size_t i = 0; i + 3 < num_runs; i += 4) {
		__m128 v = _mm_loadu_ps((float *)&x[i]);
		func(v);
	}
	uint64_t toc = __rdtsc();
	return (toc - tic);
}

uint64_t time_trig2(float (*func)(float, float), const float *x, const float *y,
		    const uint64_t num_runs) {
	/* Warm the function up */
	for (int i = 0; i < 1000; ++i) {
		func(x[i], y[i]);
	}

	uint64_t tic = __rdtsc();
	for (uint64_t i = 0; i < num_runs; ++i) {

		func(x[i], y[i]);
	}
	uint64_t toc = __rdtsc();

	return (toc - tic);
}

uint64_t time_trig2_sse(__m128 (*func)(__m128, __m128), const float *x,
			const float *y, const uint64_t num_runs) {
	uint64_t tic = __rdtsc();

	for (size_t i = 0; i + 3 < num_runs; i += 4) {
		__m128 v = _mm_loadu_ps((float *)&x[i]);
		__m128 w = _mm_loadu_ps((float *)&y[i]);
		func(v, w);
	}
	uint64_t toc = __rdtsc();
	return (toc - tic);
}

bool float_within(float f1, float f2, float delta) {
	return fabs(f1 - f2) <= delta;
}

static __m128 mask_to_neg1_or_1(__m128 mask) {
	__m128 ones = _mm_set1_ps(1.0f);
	__m128 negs = _mm_set1_ps(-2.0f);
	__m128 flip = _mm_and_ps(mask, negs);
	return _mm_add_ps(ones, flip);
}

__m128 sin_sse(__m128 x) {
	const __m128 inv_2pi	  = _mm_set1_ps(_MM_1_OVER_2PI);
	const __m128 two_pi	  = _mm_set1_ps(_MM_2PI);
	const __m128 pi		  = _mm_set1_ps(_MM_PI);
	const __m128 pi_2	  = _mm_set1_ps(_MM_PI_2);
	const __m128 pi_2_times_3 = _mm_set1_ps(_MM_PI_2 * 3.0f);

	__m128 k     = _mm_mul_ps(x, inv_2pi);
	k	     = _mm_round_ps(k, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
	__m128 x_mod = _mm_sub_ps(x, _mm_mul_ps(k, two_pi));

	__m128 in_q1 = _mm_cmple_ps(x_mod, pi_2);
	__m128 in_q2 =
	    _mm_and_ps(_mm_cmpgt_ps(x_mod, pi_2), _mm_cmple_ps(x_mod, pi));
	__m128 in_q3 = _mm_and_ps(_mm_cmpgt_ps(x_mod, pi),
				  _mm_cmple_ps(x_mod, pi_2_times_3));
	__m128 in_q4 = _mm_cmpgt_ps(x_mod, pi_2_times_3);

	__m128 x_1	 = _mm_and_ps(in_q1, x_mod);
	__m128 x_2	 = _mm_and_ps(in_q2, _mm_sub_ps(pi, x_mod));
	__m128 x_3	 = _mm_and_ps(in_q3, _mm_sub_ps(x_mod, pi));
	__m128 x_4	 = _mm_and_ps(in_q4, _mm_sub_ps(two_pi, x_mod));
	__m128 x_reduced = _mm_or_ps(_mm_or_ps(x_1, x_2), _mm_or_ps(x_3, x_4));

	__m128 flip = _mm_or_ps(in_q3, in_q4);
	__m128 sign = mask_to_neg1_or_1(flip);

	__m128 x2 = _mm_mul_ps(x_reduced, x_reduced);
	__m128 y  = _mm_set1_ps(-0.00000002f);
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(0.00000275f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(-0.00019841f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(0.00833333f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(-0.16666667f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(1.00000000f));

	y = _mm_mul_ps(y, x_reduced);
	return _mm_mul_ps(y, sign);
}

__m128 cos_sse(__m128 x) {
	const __m128 inv_2pi	  = _mm_set1_ps(_MM_1_OVER_2PI);
	const __m128 two_pi	  = _mm_set1_ps(_MM_2PI);
	const __m128 pi		  = _mm_set1_ps(_MM_PI);
	const __m128 pi_2	  = _mm_set1_ps(_MM_PI_2);
	const __m128 pi_2_times_3 = _mm_set1_ps(_MM_PI_2 * 3.0f);

	__m128 k     = _mm_mul_ps(x, inv_2pi);
	k	     = _mm_round_ps(k, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
	__m128 x_mod = _mm_sub_ps(x, _mm_mul_ps(k, two_pi));

	__m128 in_q1 = _mm_cmple_ps(x_mod, pi_2);
	__m128 in_q2 =
	    _mm_and_ps(_mm_cmpgt_ps(x_mod, pi_2), _mm_cmple_ps(x_mod, pi));
	__m128 in_q3 = _mm_and_ps(_mm_cmpgt_ps(x_mod, pi),
				  _mm_cmple_ps(x_mod, pi_2_times_3));
	__m128 in_q4 = _mm_cmpgt_ps(x_mod, pi_2_times_3);

	__m128 x1	 = _mm_and_ps(in_q1, x_mod);
	__m128 x2	 = _mm_and_ps(in_q2, _mm_sub_ps(pi, x_mod));
	__m128 x3	 = _mm_and_ps(in_q3, _mm_sub_ps(x_mod, pi));
	__m128 x4	 = _mm_and_ps(in_q4, _mm_sub_ps(two_pi, x_mod));
	__m128 x_reduced = _mm_or_ps(_mm_or_ps(x1, x2), _mm_or_ps(x3, x4));

	__m128 flip = _mm_or_ps(in_q2, in_q3);
	__m128 sign = mask_to_neg1_or_1(flip);

	__m128 x2_val = _mm_mul_ps(x_reduced, x_reduced);
	__m128 y      = _mm_set1_ps(-0.000000263435f);
	y	      = _mm_fmadd_ps(y, x2_val, _mm_set1_ps(0.000024775785f));
	y	      = _mm_fmadd_ps(y, x2_val, _mm_set1_ps(-0.001388864030f));
	y	      = _mm_fmadd_ps(y, x2_val, _mm_set1_ps(0.041666655661f));
	y	      = _mm_fmadd_ps(y, x2_val, _mm_set1_ps(-0.499999998076f));
	y	      = _mm_fmadd_ps(y, x2_val, _mm_set1_ps(0.999999999921f));

	return _mm_mul_ps(y, sign);
}

__m128 atan_sse(__m128 x) {
	const __m128 one  = _mm_set1_ps(1.0f);
	const __m128 pi_2 = _mm_set1_ps(_MM_PI_2);

	// Absolute value and sign mask
	__m128 sign_mask =
	    _mm_and_ps(x, _mm_set1_ps(-0.0f));	    // Keep only sign bit
	__m128 abs_x = _mm_andnot_ps(sign_mask, x); // Clear sign bit

	// Range reduction mask: if x > 1, use identity atan(x) = π/2 −
	// atan(1/x)
	__m128 use_identity_mask = _mm_cmpgt_ps(abs_x, one);

	// Compute reciprocal for inputs where |x| > 1
	__m128 x_inv = _mm_div_ps(one, abs_x);

	// Select reduced x: either x or 1/x
	__m128 x_red = _mm_or_ps(_mm_and_ps(use_identity_mask, x_inv),
				 _mm_andnot_ps(use_identity_mask, abs_x));

	// Polynomial eval on x_red²
	__m128 x2 = _mm_mul_ps(x_red, x_red);
	__m128 y  = _mm_set1_ps(-0.01891985f);
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(0.06909289f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(-0.12945813f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(0.19787976f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(-0.33319414f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(0.99999750f));

	y = _mm_mul_ps(y, x_red); // atan is odd

	// If we used identity, atan(x) = π/2 − atan(1/x)
	__m128 atan_corrected = _mm_sub_ps(pi_2, y);
	y = _mm_or_ps(_mm_and_ps(use_identity_mask, atan_corrected),
		      _mm_andnot_ps(use_identity_mask, y));

	// Restore original sign
	y = _mm_or_ps(y, sign_mask);
	return y;
}

__m128 atan2_sse(__m128 y, __m128 x) {
	const __m128 zero     = _mm_set1_ps(0.0f);
	const __m128 pi	      = _mm_set1_ps(_MM_PI);
	const __m128 pi_2     = _mm_set1_ps(_MM_PI_2);
	const __m128 neg_pi_2 = _mm_set1_ps(-_MM_PI_2);

	// Compute atan(y/x) with correct handling of division by zero
	__m128 abs_x = _mm_andnot_ps(_mm_set1_ps(-0.0f), x); // fabs(x)
	__m128 abs_y = _mm_andnot_ps(_mm_set1_ps(-0.0f), y); // fabs(y)

	// Prevent divide by zero by setting zero denom to 1 (won't matter due
	// to masking)
	__m128 safe_x =
	    _mm_blendv_ps(x, _mm_set1_ps(1.0f), _mm_cmpeq_ps(abs_x, zero));
	__m128 z      = _mm_div_ps(y, safe_x);
	__m128 atan_z = atan_sse(z);

	// Masks for quadrant handling
	__m128 x_gt0 = _mm_cmpgt_ps(x, zero); // x > 0
	__m128 x_lt0 = _mm_cmplt_ps(x, zero); // x < 0
	__m128 y_eq0 = _mm_cmpeq_ps(y, zero); // y == 0
	__m128 y_lt0 = _mm_cmplt_ps(y, zero); // y < 0
	__m128 y_gt0 = _mm_cmpgt_ps(y, zero); // y > 0
	__m128 x_eq0 = _mm_cmpeq_ps(x, zero); // x == 0

	// atan2(y,x) = atan(y/x) if x > 0
	//            = atan(y/x) + π if x < 0 and y >= 0
	//            = atan(y/x) - π if x < 0 and y < 0
	__m128 add_pi = _mm_and_ps(x_lt0, y_gt0);
	__m128 sub_pi = _mm_and_ps(x_lt0, y_lt0);
	__m128 result = atan_z;
	result	      = _mm_add_ps(result, _mm_and_ps(add_pi, pi));
	result	      = _mm_sub_ps(result, _mm_and_ps(sub_pi, pi));

	// Handle x == 0
	__m128 x0_y0   = _mm_and_ps(x_eq0, y_eq0); // atan2(0,0) → 0
	__m128 x0_ygt0 = _mm_and_ps(x_eq0, y_gt0); // atan2(+y,0) → +π/2
	__m128 x0_ylt0 = _mm_and_ps(x_eq0, y_lt0); // atan2(-y,0) → -π/2

	result = _mm_blendv_ps(result, zero, x0_y0);
	result = _mm_blendv_ps(result, pi_2, x0_ygt0);
	result = _mm_blendv_ps(result, neg_pi_2, x0_ylt0);

	return result;
}

__m128 sinlr_sse(__m128 x) {
	const __m128 inv_2pi = _mm_set1_ps(_MM_1_OVER_2PI);
	const __m128 two_pi  = _mm_set1_ps(_MM_2PI);
	const __m128 pi	     = _mm_set1_ps(_MM_PI);

	__m128 k = _mm_mul_ps(x, inv_2pi);
	k	 = _mm_round_ps(k, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
	x	 = _mm_sub_ps(x, _mm_mul_ps(k, two_pi));
	x	 = _mm_sub_ps(x, pi);

	__m128 x2 = _mm_mul_ps(x, x);
	__m128 y  = _mm_set1_ps(0.00000000f);
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(-0.00000002f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(0.00000275f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(-0.00019841f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(0.00833333f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(-0.16666666f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(1.00000000f));

	y = _mm_mul_ps(y, -x);
	return _mm_min_ps(_mm_max_ps(y, _mm_set1_ps(-1.0f)), _mm_set1_ps(1.0f));
}

__m128 coslr_sse(__m128 x) {
	const __m128 inv_2pi = _mm_set1_ps(_MM_1_OVER_2PI);
	const __m128 two_pi  = _mm_set1_ps(_MM_2PI);
	const __m128 pi	     = _mm_set1_ps(_MM_PI);

	__m128 k = _mm_mul_ps(x, inv_2pi);
	k	 = _mm_round_ps(k, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
	x	 = _mm_sub_ps(x, _mm_mul_ps(k, two_pi));
	x	 = _mm_sub_ps(x, pi);

	__m128 x2 = _mm_mul_ps(x, x);
	__m128 y  = _mm_set1_ps(0.00000000f);
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(-0.00000027f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(0.00002478f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(-0.00138884f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(0.04166661f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(-0.49999997f));
	y	  = _mm_fmadd_ps(y, x2, _mm_set1_ps(1.00000000f));

	return _mm_mul_ps(y, _mm_set1_ps(-1.0f));
}

int main() {
	uint64_t not_within = 0;
	srand(__rdtsc());
	float *x = malloc(sizeof(float) * RUNS);
	for (size_t i = 0; i < RUNS; ++i) {
	  x[i] = (float)rand() / (float)rand() / (2.0f * M_PI * 16.0f);
	}

	uint64_t real_sin_time = time_trig(sinf, x, RUNS);
	uint64_t sin_sse_time  = time_trig_sse(sin_sse, x, RUNS);

	for (size_t i = 0; i + 3 < RUNS; i += 4) {
		__m128 v	= _mm_loadu_ps((float *)&x[i]);
		__m128 sse_sins = sin_sse(v);
		float  sins[4];
		_mm_storeu_ps(sins, sse_sins);
		for (size_t s = 0; s < 4; ++s) {
			if (!float_within(sinf(x[i + s]), sins[s], DELTA)) {
				fprintf(stderr,
					"x: %fpi, sin: %f, sin_sse: %f\n",
					fmod(x[i + s] / M_PI, 2.0f * M_PI),
					sinf(x[i + s]), sins[s]);
				++not_within;
			}
		}
	}
	printf("sin_sse takes %f times as long as sinf\n",
	       (double)sin_sse_time / (double)real_sin_time);
	printf("sin_sse was within DELTA %f%% of the time.\n",
	       ((float)(RUNS - not_within) / (float)RUNS) * 100.0f);

	/***********************************************************************
	 **********************************************************************/

	not_within = 0.0f;
	for (size_t i = 0; i < RUNS; ++i) {
	  x[i] = (float)rand() / (float)rand() * (2.0f * M_PI * 16.0f);
	}

	uint64_t real_cos_time = time_trig(cosf, x, RUNS);
	uint64_t cos_sse_time  = time_trig_sse(cos_sse, x, RUNS);

	for (size_t i = 0; i + 3 < RUNS; i += 4) {
		__m128 v	= _mm_loadu_ps((float *)&x[i]);
		__m128 sse_coss = cos_sse(v);
		float  coss[4];
		_mm_storeu_ps(coss, sse_coss);
		for (size_t s = 0; s < 4; ++s) {
			if (!float_within(cosf(x[i + s]), coss[s], DELTA)) {
				fprintf(stderr,
					"x: %fpi, cos: %f, cos_sse: %f\n",
					fmod(x[i + s] / M_PI, 2.0f * M_PI),
					cosf(x[i + s]), coss[s]);
				++not_within;
			}
		}
	}
	printf("cos_sse takes %f times as long as cosf\n",
	       (double)cos_sse_time / (double)real_cos_time);
	printf("cos_sse was within DELTA %f%% of the time.\n",
	       ((float)(RUNS - not_within) / (float)RUNS) * 100.0f);

	/***********************************************************************
	 **********************************************************************/

	not_within = 0.0f;
	float *y   = malloc(sizeof(float) * RUNS);
	for (size_t i = 0; i < RUNS; ++i) {
	  y[i] = (float)rand() / (float)rand() / (2 * M_PI * 16.0f);
	}

	uint64_t real_atan2_time = time_trig2(atan2f, x, y, RUNS);
	uint64_t atan2_sse_time	 = time_trig2_sse(atan2_sse, x, y, RUNS);

	for (size_t i = 0; i + 3 < RUNS; i += 4) {
		__m128 v	  = _mm_loadu_ps((float *)&x[i]);
		__m128 w	  = _mm_loadu_ps((float *)&y[i]);
		__m128 sse_atan2s = atan2_sse(v, w);
		float  atan2s[4];
		_mm_storeu_ps(atan2s, sse_atan2s);
		for (size_t s = 0; s < 4; ++s) {
			if (!float_within(atan2f(x[i + s], y[i + s]), atan2s[s],
					  DELTA)) {
				fprintf(stderr,
					"x: %fpi, atan2: %f, atan2_sse: %f\n",
					fmod(x[i + s] / M_PI, 2.0f * M_PI),
					atan2f(x[i + s], y[i + s]), atan2s[s]);
				++not_within;
			}
		}
	}
	printf("atan2_sse takes %f times as long as atan2f\n",
	       (double)atan2_sse_time / (double)real_atan2_time);
	printf("atan2_sse was within DELTA %f%% of the time.\n",
	       ((float)(RUNS - not_within) / (float)RUNS) * 100.0f);

	/***********************************************************************
	 **********************************************************************/

	not_within = 0.0f;
	for (size_t i = 0; i < RUNS; ++i) {
	  x[i] = (float)rand() / (float)rand() / (2 * M_PI * 16.0f);
	}

	uint64_t real_sinlr_time = time_trig(sinf, x, RUNS);
	uint64_t sinlr_sse_time	 = time_trig_sse(sinlr_sse, x, RUNS);

	for (size_t i = 0; i + 3 < RUNS; i += 4) {
		__m128 v	  = _mm_loadu_ps((float *)&x[i]);
		__m128 sse_sinlrs = sinlr_sse(v);
		float  sinlrs[4];
		_mm_storeu_ps(sinlrs, sse_sinlrs);
		for (size_t s = 0; s < 4; ++s) {
			if (!float_within(sinf(x[i + s]), sinlrs[s], DELTA)) {
				fprintf(stderr,
					"x: %f, reduced: %fpi, sin: %f, "
					"sinlr_sse: %f\n",
					x[i + s],
					fmod(x[i + s], 2.0f * M_PI) / M_PI,
					sinf(x[i + s]), sinlrs[s]);
				++not_within;
			}
		}
	}
	printf("sinlr_sse takes %f times as long as sinf\n",
	       (double)sinlr_sse_time / (double)real_sinlr_time);
	printf("sinlr_sse was within DELTA %f%% of the time.\n",
	       ((float)(RUNS - not_within) / (float)RUNS) * 100.0f);

	/***********************************************************************
	 **********************************************************************/

	not_within = 0.0f;
	for (size_t i = 0; i < RUNS; ++i) {
	  x[i] = (float)rand() / (float)rand() / (2 * M_PI * 16.0f);
	}

	uint64_t real_coslr_time = time_trig(cosf, x, RUNS);
	uint64_t coslr_sse_time	 = time_trig_sse(coslr_sse, x, RUNS);

	for (size_t i = 0; i + 3 < RUNS; i += 4) {
		__m128 v	  = _mm_loadu_ps((float *)&x[i]);
		__m128 sse_coslrs = coslr_sse(v);
		float  coslrs[4];
		_mm_storeu_ps(coslrs, sse_coslrs);
		for (size_t s = 0; s < 4; ++s) {
			if (!float_within(cosf(x[i + s]), coslrs[s], DELTA)) {
				fprintf(stderr,
					"x: %f, reduced: %fpi, cos: %f, "
					"coslr_sse: %f\n",
					x[i + s],
					fmod(x[i + s], 2.0f * M_PI) / M_PI,
					cosf(x[i + s]), coslrs[s]);
				++not_within;
			}
		}
	}
	printf("coslr_sse takes %f times as long as cosf\n",
	       (double)coslr_sse_time / (double)real_coslr_time);
	printf("coslr_sse was within DELTA %f%% of the time.\n",
	       ((float)(RUNS - not_within) / (float)RUNS) * 100.0f);

	/***********************************************************************
	 **********************************************************************/

	free(x);
	free(y);

	return 0;
}
