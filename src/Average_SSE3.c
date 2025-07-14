#ifdef PALETTE_DEBUG
#include <assert.h>
#endif
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#if defined(__x86_64__) || defined(_M_X64)
#include <pmmintrin.h>
#endif

#include "Color.h"

static inline float horizontal_sum(__m128 v) {
	__m128 shuf = _mm_movehdup_ps(v);
	__m128 sums = _mm_add_ps(v, shuf);
	shuf	    = _mm_movehl_ps(shuf, sums);
	sums	    = _mm_add_ss(sums, shuf);
	return _mm_cvtss_f32(sums);
}

struct Color cielab_avg_sse3(const struct cielab_SoA *colors,
			     uint16_t		      num_col) {
	__m128 sum_l = _mm_setzero_ps();
	__m128 sum_a = _mm_setzero_ps();
	__m128 sum_b = _mm_setzero_ps();

	size_t i = 0;
	for (; i + 15 < num_col; i += 16) {
#ifdef PALETTE_DEBUG
		assert(((uintptr_t)&colors->l[i + 0] % 16) == 0);
		assert(((uintptr_t)&colors->a[i + 0] % 16) == 0);
		assert(((uintptr_t)&colors->b[i + 0] % 16) == 0);

		assert(((uintptr_t)&colors->l[i + 4] % 16) == 0);
		assert(((uintptr_t)&colors->a[i + 4] % 16) == 0);
		assert(((uintptr_t)&colors->b[i + 4] % 16) == 0);

		assert(((uintptr_t)&colors->l[i + 8] % 16) == 0);
		assert(((uintptr_t)&colors->a[i + 8] % 16) == 0);
		assert(((uintptr_t)&colors->b[i + 8] % 16) == 0);

		assert(((uintptr_t)&colors->l[i + 12] % 16) == 0);
		assert(((uintptr_t)&colors->a[i + 12] % 16) == 0);
		assert(((uintptr_t)&colors->b[i + 12] % 16) == 0);
#endif
		__m128 l1 = _mm_load_ps(&colors->l[i + 0]);
		__m128 l2 = _mm_load_ps(&colors->l[i + 4]);
		__m128 l3 = _mm_load_ps(&colors->l[i + 8]);
		__m128 l4 = _mm_load_ps(&colors->l[i + 12]);

		__m128 a1 = _mm_load_ps(&colors->a[i + 0]);
		__m128 a2 = _mm_load_ps(&colors->a[i + 4]);
		__m128 a3 = _mm_load_ps(&colors->a[i + 8]);
		__m128 a4 = _mm_load_ps(&colors->a[i + 12]);

		__m128 b1 = _mm_load_ps(&colors->b[i + 0]);
		__m128 b2 = _mm_load_ps(&colors->b[i + 4]);
		__m128 b3 = _mm_load_ps(&colors->b[i + 8]);
		__m128 b4 = _mm_load_ps(&colors->b[i + 12]);

		sum_l = _mm_add_ps(sum_l, l1);
		sum_l = _mm_add_ps(sum_l, l2);
		sum_l = _mm_add_ps(sum_l, l3);
		sum_l = _mm_add_ps(sum_l, l4);

		sum_a = _mm_add_ps(sum_a, a1);
		sum_a = _mm_add_ps(sum_a, a2);
		sum_a = _mm_add_ps(sum_a, a3);
		sum_a = _mm_add_ps(sum_a, a4);

		sum_b = _mm_add_ps(sum_b, b1);
		sum_b = _mm_add_ps(sum_b, b2);
		sum_b = _mm_add_ps(sum_b, b3);
		sum_b = _mm_add_ps(sum_b, b4);
	}
	float l, a, b;
	l = a = b = 0.0f;
	float l_vec[4];
	_mm_storeu_ps(l_vec, sum_l);
	l = l_vec[0] + l_vec[1] + l_vec[2] + l_vec[3];
	float a_vec[4];
	_mm_storeu_ps(a_vec, sum_a);
	l = a_vec[0] + a_vec[1] + a_vec[2] + a_vec[3];
	float b_vec[4];
	_mm_storeu_ps(b_vec, sum_b);
	l = b_vec[0] + b_vec[1] + b_vec[2] + b_vec[3];

	for (; i < num_col; ++i) {
		l += colors->l[i];
		a += colors->a[i];
		b += colors->b[i];
	}

	float scale = 1.0f / num_col;

	return Color_create(l * scale, a * scale, b * scale);
}
