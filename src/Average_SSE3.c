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

struct Color cielab_avg_sse3(const struct cielab_SoA *__restrict colors,
			     uint16_t		      num_col) {
	__m128 sum_l1 = _mm_setzero_ps();
	__m128 sum_a1 = _mm_setzero_ps();
	__m128 sum_b1 = _mm_setzero_ps();

	__m128 sum_l2 = _mm_setzero_ps();
	__m128 sum_a2 = _mm_setzero_ps();
	__m128 sum_b2 = _mm_setzero_ps();

	__m128 sum_l3 = _mm_setzero_ps();
	__m128 sum_a3 = _mm_setzero_ps();
	__m128 sum_b3 = _mm_setzero_ps();

	__m128 sum_l4 = _mm_setzero_ps();
	__m128 sum_a4 = _mm_setzero_ps();
	__m128 sum_b4 = _mm_setzero_ps();

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

		sum_l1 = _mm_add_ps(sum_l1, l1);
		sum_l2 = _mm_add_ps(sum_l2, l2);
		sum_l3 = _mm_add_ps(sum_l3, l3);
		sum_l4 = _mm_add_ps(sum_l4, l4);

		sum_a1 = _mm_add_ps(sum_a1, a1);
		sum_a2 = _mm_add_ps(sum_a2, a2);
		sum_a3 = _mm_add_ps(sum_a3, a3);
		sum_a4 = _mm_add_ps(sum_a4, a4);

		sum_b1 = _mm_add_ps(sum_b1, b1);
		sum_b2 = _mm_add_ps(sum_b2, b2);
		sum_b3 = _mm_add_ps(sum_b3, b3);
		sum_b4 = _mm_add_ps(sum_b4, b4);
	}

	__m128 sum_l =
	    _mm_add_ps(_mm_add_ps(sum_l1, sum_l2), _mm_add_ps(sum_l3, sum_l4));
	__m128 sum_a =
	    _mm_add_ps(_mm_add_ps(sum_a1, sum_a2), _mm_add_ps(sum_a3, sum_a4));
	__m128 sum_b =
	    _mm_add_ps(_mm_add_ps(sum_b1, sum_b2), _mm_add_ps(sum_b3, sum_b4));

	float l_vec[4];
	_mm_store_ps(l_vec, sum_l);
	float l = l_vec[0] + l_vec[1] + l_vec[2] + l_vec[3];

	float a_vec[4];
	_mm_store_ps(a_vec, sum_a);
	float a = a_vec[0] + a_vec[1] + a_vec[2] + a_vec[3];

	float b_vec[4];
	_mm_store_ps(b_vec, sum_b);
	float b = b_vec[0] + b_vec[1] + b_vec[2] + b_vec[3];

	for (; i < num_col; ++i) {
		l += colors->l[i];
		a += colors->a[i];
		b += colors->b[i];
	}

	float scale = 1.0f / num_col;

	return Color_create_cielab(l * scale, a * scale, b * scale);
}

struct Color oklab_avg_sse3(const struct oklab_SoA *__restrict colors,
			     uint16_t		      num_col) {
	__m128 sum_l1 = _mm_setzero_ps();
	__m128 sum_a1 = _mm_setzero_ps();
	__m128 sum_b1 = _mm_setzero_ps();

	__m128 sum_l2 = _mm_setzero_ps();
	__m128 sum_a2 = _mm_setzero_ps();
	__m128 sum_b2 = _mm_setzero_ps();

	__m128 sum_l3 = _mm_setzero_ps();
	__m128 sum_a3 = _mm_setzero_ps();
	__m128 sum_b3 = _mm_setzero_ps();

	__m128 sum_l4 = _mm_setzero_ps();
	__m128 sum_a4 = _mm_setzero_ps();
	__m128 sum_b4 = _mm_setzero_ps();

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

		sum_l1 = _mm_add_ps(sum_l1, l1);
		sum_l2 = _mm_add_ps(sum_l2, l2);
		sum_l3 = _mm_add_ps(sum_l3, l3);
		sum_l4 = _mm_add_ps(sum_l4, l4);

		sum_a1 = _mm_add_ps(sum_a1, a1);
		sum_a2 = _mm_add_ps(sum_a2, a2);
		sum_a3 = _mm_add_ps(sum_a3, a3);
		sum_a4 = _mm_add_ps(sum_a4, a4);

		sum_b1 = _mm_add_ps(sum_b1, b1);
		sum_b2 = _mm_add_ps(sum_b2, b2);
		sum_b3 = _mm_add_ps(sum_b3, b3);
		sum_b4 = _mm_add_ps(sum_b4, b4);
	}

	__m128 sum_l =
	    _mm_add_ps(_mm_add_ps(sum_l1, sum_l2), _mm_add_ps(sum_l3, sum_l4));
	__m128 sum_a =
	    _mm_add_ps(_mm_add_ps(sum_a1, sum_a2), _mm_add_ps(sum_a3, sum_a4));
	__m128 sum_b =
	    _mm_add_ps(_mm_add_ps(sum_b1, sum_b2), _mm_add_ps(sum_b3, sum_b4));

	float l_vec[4];
	_mm_store_ps(l_vec, sum_l);
	float l = l_vec[0] + l_vec[1] + l_vec[2] + l_vec[3];

	float a_vec[4];
	_mm_store_ps(a_vec, sum_a);
	float a = a_vec[0] + a_vec[1] + a_vec[2] + a_vec[3];

	float b_vec[4];
	_mm_store_ps(b_vec, sum_b);
	float b = b_vec[0] + b_vec[1] + b_vec[2] + b_vec[3];

	for (; i < num_col; ++i) {
		l += colors->l[i];
		a += colors->a[i];
		b += colors->b[i];
	}

	float scale = 1.0f / num_col;

	return Color_create_oklab(l * scale, a * scale, b * scale);
}
