#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <immintrin.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <x86intrin.h>

#define DELTA 1e-1

#define FMADD(a, b, c) _mm_add_ps(_mm_mul_ps((a), (b)), (c))
#define _mm_abs_ps(x) _mm_and_ps((x), _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF)))
#define _MM_PI 3.14159265f
#define _MM_2PI 6.28318531f
#define _MM_PI_2 1.57079633f
#define _MM_1_OVER_2PI 0.159154943f
#define FL_PER_VEC 4

static const int RUNS = 65536;

static inline float clamp(float x, float min, float max) {
	return x < min ? min : (x > max ? max : x);
}

static inline __m128 V(float x) {
	return _mm_set1_ps(x);
}

static inline int time_scalar_math(float (*func)(float), float *nums) {
	volatile float total = 0;
	/* Warm the function up */
	for (int i = 0; i < 10000; ++i) {
		total += func(nums[i]);
	}

	uint64_t tic = __rdtsc();
	for (int i = 0; i < RUNS; ++i) {

		total += func(nums[i]);
	}
	uint64_t toc = __rdtsc();

	return (toc - tic) / RUNS;
}

static inline int time_scalar_math2(float (*func)(float, float), float *nums1,
				    float *nums2) {
	volatile float total = 0;
	/* Warm the function up */
	for (int i = 0; i < 10000; ++i) {
		total += func(nums1[i], nums2[i]);
	}

	uint64_t tic = __rdtsc();
	for (int i = 0; i < RUNS; ++i) {

		total += func(nums1[i], nums2[i]);
	}
	uint64_t toc = __rdtsc();

	return (toc - tic) / RUNS;
}

static inline int time_vector_math(__m128 (*func)(__m128), __m128 *nums) {
	volatile __m128 total = _mm_setzero_ps();
	/* Warm up the function */
	for (int i = 0; i < 10000; ++i) {
		func(nums[i]);
	}

	uint64_t tic = __rdtsc();
	for (int i = 0; i < RUNS / FL_PER_VEC; ++i) {
		total = _mm_add_ps(total, func(nums[i]));
	}
	uint64_t toc = __rdtsc();

	return (toc - tic) / RUNS;
}

static inline int time_vector_math2(__m128 (*func)(__m128, __m128),
				    __m128 *nums1, __m128 *nums2) {
	volatile __m128 total = _mm_setzero_ps();
	/* Warm up the function */
	for (int i = 0; i < 10000; ++i) {
		func(nums1[i], nums2[i]);
	}

	uint64_t tic = __rdtsc();
	for (int i = 0; i < RUNS / FL_PER_VEC; ++i) {
		total = _mm_add_ps(total, func(nums1[i], nums2[i]));
	}
	uint64_t toc = __rdtsc();

	return (toc - tic) / RUNS;
}

static inline float tmaxf(float a, float b) {
	return (a + b + fabs(a - b)) / 2;
}

static inline float test_math_accuracy(float (*sc_fn)(float), float    *nums,
				       __m128 (*ve_fn)(__m128), __m128 *vecs) {
	float max_err	     = 0.0f;
	float vec_answers[4] = {0};

	for (int j = 0, i = 0; j + FL_PER_VEC <= RUNS; j += 4, ++i) {
		__m128 cur_ans = ve_fn(vecs[i]);
		_mm_storeu_ps(vec_answers, cur_ans);
		for (int k = 0; k < FL_PER_VEC; ++k) {
			max_err = tmaxf(max_err, fabsf(vec_answers[k] -
						       sc_fn(nums[j + k])));
		}
	}
	return max_err;
}

static inline float test_math_accuracy2(float (*sc_fn)(float, float),
					float *nums1, float *nums2,
					__m128 (*ve_fn)(__m128, __m128),
					__m128 *vecs1, __m128 *vecs2) {
	float max_err	     = 0.0f;
	float vec_answers[4] = {0};

	for (int j = 0, i = 0; j + FL_PER_VEC <= RUNS; j += 4, ++i) {
		__m128 cur_ans = ve_fn(vecs1[i], vecs2[i]);
		_mm_storeu_ps(vec_answers, cur_ans);
		for (int k = 0; k < FL_PER_VEC; ++k) {
			max_err = tmaxf(
			    max_err, fabsf(vec_answers[k] -
					   sc_fn(nums1[j + k], nums2[j + k])));
		}
	}
	return max_err;
}

static inline bool float_within(float f1, float f2, float delta) {
	return fabs(f1 - f2) <= delta;
}

static inline __m128 _mm_atan_ps(__m128 x) {
	const __m128 one  = V(1.0f);
	const __m128 pi_2 = V(_MM_PI_2);

	__m128 sign_mask = _mm_and_ps(x, V(-0.0f));
	__m128 abs_x	 = _mm_andnot_ps(sign_mask, x);

	__m128 use_identity_mask = _mm_cmpgt_ps(abs_x, one);

	__m128 x_inv = _mm_div_ps(one, abs_x);

	x = _mm_or_ps(_mm_and_ps(use_identity_mask, x_inv),
		      _mm_andnot_ps(use_identity_mask, abs_x));

	
	__m128 coeffs[6] = {
	    _mm_set1_ps(-1.889933876e-3f), _mm_set1_ps(1.723863386e-2f),
	    _mm_set1_ps(-6.777403981e-2f), _mm_set1_ps(0.163024798f),
	    _mm_set1_ps(-0.324574113f),	   _mm_set1_ps(0.999377847f)};
	
	__m128 y  = coeffs[0];
	__m128 x2 = _mm_mul_ps(x, x);
	y	  = FMADD(y, x2, coeffs[1]);
	y	  = FMADD(y, x2, coeffs[2]);
	y	  = FMADD(y, x2, coeffs[3]);
	y	  = FMADD(y, x2, coeffs[4]);
	y	  = FMADD(y, x2, coeffs[5]);

	y = _mm_mul_ps(y, x);

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

	/* Range reduce to -pi/2, pi/2 */
	__m128 k = _mm_mul_ps(x, inv_2pi);
	k	 = _mm_round_ps(k, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
	x	 = _mm_sub_ps(x, _mm_mul_ps(k, two_pi));
	x	 = _mm_sub_ps(x, pi);

	__m128 coeffs[6] = {
	    _mm_set1_ps(-2.379471354e-8f), _mm_set1_ps(2.751885564e-6f),
	    _mm_set1_ps(-1.984070286e-4f), _mm_set1_ps(8.333329264e-3f),
	    _mm_set1_ps(-0.166666672f),	   _mm_set1_ps(0.999999940f)};

	__m128 x2 = _mm_mul_ps(x, x);
	__m128 y  = coeffs[0];

	y = FMADD(y, x2, coeffs[1]);
	y = FMADD(y, x2, coeffs[2]);
	y = FMADD(y, x2, coeffs[3]);
	y = FMADD(y, x2, coeffs[4]);
	y = FMADD(y, x2, coeffs[5]);

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

	__m128 coeffs[6] = {
	    _mm_set1_ps(-2.605149520e-7f),  _mm_set1_ps(2.476016135e-5f),
	    _mm_set1_ps(-1.388836140e-03f), _mm_set1_ps(4.166663625e-02f),
	    _mm_set1_ps(-0.499999993f),	    _mm_set1_ps(0.999999999)};
	__m128 x2 = _mm_mul_ps(x, x);
	__m128 y  = coeffs[0];
	y	  = FMADD(y, x2, coeffs[1]);
	y	  = FMADD(y, x2, coeffs[2]);
	y	  = FMADD(y, x2, coeffs[3]);
	y	  = FMADD(y, x2, coeffs[4]);
	y	  = FMADD(y, x2, coeffs[5]);

	y = _mm_mul_ps(y, V(-1.0f));
	return _mm_min_ps(_mm_max_ps(y, V(-1.0f)), V(1.0f));
}

static inline __m128 _mm_cbrt_ps(__m128 x) {
	__m128 coeffs[8] = {
	    _mm_set1_ps(527.924255f),  _mm_set1_ps(-1886.101562f),
	    _mm_set1_ps(2668.042236f), _mm_set1_ps(-1896.805298f),
	    _mm_set1_ps(710.645935f),  _mm_set1_ps(-135.148804f),
	    _mm_set1_ps(12.443144f),   _mm_set1_ps(5.604970038e-2f)};

	__m128 y = coeffs[0];
	y	 = FMADD(y, x, coeffs[1]);
	y	 = FMADD(y, x, coeffs[2]);
	y	 = FMADD(y, x, coeffs[3]);
	y	 = FMADD(y, x, coeffs[4]);
	y	 = FMADD(y, x, coeffs[5]);
	y	 = FMADD(y, x, coeffs[6]);
	y	 = FMADD(y, x, coeffs[7]);

	return y;
}

static inline __m128 _mm_exp_ps(__m128 x) {
	__m128 coeffs[9] = {
	    _mm_set1_ps(1.956499731e-8f), _mm_set1_ps(1.677241844e-6f),
	    _mm_set1_ps(6.017330956e-5f), _mm_set1_ps(1.174603360e-3f),
	    _mm_set1_ps(1.359774692e-2f), _mm_set1_ps(9.559546303e-2f),
	    _mm_set1_ps(0.401840591f),	  _mm_set1_ps(0.945352098f),
	    _mm_set1_ps(0.994917907f),
	};

	__m128 y = coeffs[0];
	y	 = FMADD(y, x, coeffs[1]);
	y	 = FMADD(y, x, coeffs[2]);
	y	 = FMADD(y, x, coeffs[3]);
	y	 = FMADD(y, x, coeffs[4]);
	y	 = FMADD(y, x, coeffs[5]);
	y	 = FMADD(y, x, coeffs[6]);
	y	 = FMADD(y, x, coeffs[7]);
	y	 = FMADD(y, x, coeffs[8]);

	return y;
}

int main() {
	srand(__rdtsc());
	float  *nums  = malloc(sizeof(float) * RUNS);
	float  *nums2 = malloc(sizeof(float) * RUNS);
	__m128 *vecs  = malloc(sizeof(__m128) * RUNS / FL_PER_VEC);
	__m128 *vecs2 = malloc(sizeof(__m128) * RUNS / FL_PER_VEC);

	for (int i = 0; i < RUNS; ++i) {
		nums[i]	 = (float)rand() / (float)rand();
		nums2[i] = (float)rand() / (float)rand();
	}
	for (int i = 0, j = 0; i + FL_PER_VEC <= RUNS; i += FL_PER_VEC, ++j) {
		vecs[j]	 = _mm_loadu_ps(&nums[i]);
		vecs2[j] = _mm_loadu_ps(&nums2[i]);
	}
	const int   atan2_time	 = time_scalar_math2(atan2f, nums, nums2);
	const int   atan2_time_v = time_vector_math2(_mm_atan2_ps, vecs, vecs2);
	const float atan2_err =
	    test_math_accuracy2(atan2f, nums, nums2, _mm_atan2_ps, vecs, vecs2);

	for (int i = 0; i < RUNS; ++i) {
		nums[i] =
		    clamp(fmodf((float)rand() / (float)rand(), 1 + 1.2e-7),
			  -1.2e-7, 1 + 1.2e-7);
	}
	for (int i = 0, j = 0; i + FL_PER_VEC <= RUNS; i += FL_PER_VEC, ++j) {
		vecs[j] = _mm_loadu_ps(&nums[i]);
	}
	const int   cbrt_time	= time_scalar_math(cbrtf, nums);
	const int   cbrt_time_v = time_vector_math(_mm_cbrt_ps, vecs);
	const float cbrt_err =
	    test_math_accuracy(cbrtf, nums, _mm_cbrt_ps, vecs);

	for (int i = 0; i < RUNS; ++i) {
		nums[i] =
		    clamp(-fmodf((float)rand() / (float)rand(), 19.343585),
			  -19.343585, 0.0f);
	}
	for (int i = 0, j = 0; i + FL_PER_VEC <= RUNS; i += FL_PER_VEC, ++j) {
		vecs[j] = _mm_loadu_ps(&nums[i]);
	}
	const int   exp_time   = time_scalar_math(expf, nums);
	const int   exp_time_v = time_vector_math(_mm_exp_ps, vecs);
	const float exp_err = test_math_accuracy(expf, nums, _mm_exp_ps, vecs);

	for (int i = 0; i < RUNS; ++i) {
		nums[i] = clamp(fmodf((float)rand() / (float)rand(), M_PI / 2),
				-M_PI / 2, M_PI / 2);
	}
	for (int i = 0, j = 0; i + FL_PER_VEC <= RUNS; i += FL_PER_VEC, ++j) {
		vecs[j] = _mm_loadu_ps(&nums[i]);
	}
	const int   sin_time   = time_scalar_math(sinf, nums);
	const int   cos_time   = time_scalar_math(cosf, nums);
	const int   sin_time_v = time_vector_math(_mm_sin_ps, vecs);
	const int   cos_time_v = time_vector_math(_mm_cos_ps, vecs);
	const float sin_err = test_math_accuracy(sinf, nums, _mm_sin_ps, vecs);
	const float cos_err = test_math_accuracy(cosf, nums, _mm_cos_ps, vecs);

	printf("\tsin\t\t\tcos\t\t\texp\t\t\tcbrt\t\t\tatan2\n");
	printf("Scalar\t%d\t\t\t%d\t\t\t%d\t\t\t%d\t\t\t%d\n", sin_time,
	       cos_time, exp_time, cbrt_time, atan2_time);
	printf("Vector\t%d\t\t\t%d\t\t\t%d\t\t\t%d\t\t\t%d\n\n", sin_time_v,
	       cos_time_v, exp_time_v, cbrt_time_v, atan2_time_v);
	printf("Diff\t%fx\t\t%fx\t\t%fx\t\t%f\t\t%f\n",
	       (float)sin_time / (float)sin_time_v,
	       (float)cos_time / (float)cos_time_v,
	       (float)exp_time / (float)exp_time_v,
	       (float)cbrt_time / (float)cbrt_time_v,
	       (float)atan2_time / (float)atan2_time_v);
	printf("\n\nMax Err\t%f\t\t%f\t\t%f\t\t%f\t\t%f\n", sin_err, cos_err,
	       exp_err, cbrt_err, atan2_err);

	free(nums);
	free(vecs);
	free(nums2);
	free(vecs2);
	return 0;
}
