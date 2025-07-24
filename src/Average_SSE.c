#ifdef PALETTE_DEBUG
#include <assert.h>
#endif
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#if defined(__x86_64__) || defined(_M_X64)
#include <pmmintrin.h>
#endif

#include "Align.h"
#include "Color.h"

#define SQUARE(x) ((x) * (x))
void cielab_avg_sse(const struct Color *restrict colors, struct Color *restrict cents,
		    uint8_t which_cent) {
	uint16_t num_col = colors->num;
	__m128	 sum_l1	 = _mm_setzero_ps();
	__m128	 sum_a1	 = _mm_setzero_ps();
	__m128	 sum_b1	 = _mm_setzero_ps();

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
		assert(((uintptr_t)colors->alpha[i + 0] % 16) == 0);
		assert(((uintptr_t)colors->beta[i + 0] % 16) == 0);
		assert(((uintptr_t)colors->gamma[i + 0] % 16) == 0);

		assert(((uintptr_t)colors->alpha[i + 4] % 16) == 0);
		assert(((uintptr_t)colors->beta[i + 4] % 16) == 0);
		assert(((uintptr_t)colors->gamma[i + 4] % 16) == 0);

		assert(((uintptr_t)colors->alpha[i + 8] % 16) == 0);
		assert(((uintptr_t)colors->beta[i + 8] % 16) == 0);
		assert(((uintptr_t)colors->gamma[i + 8] % 16) == 0);

		assert(((uintptr_t)colors->alpha[i + 12] % 16) == 0);
		assert(((uintptr_t)colors->beta[i + 12] % 16) == 0);
		assert(((uintptr_t)colors->gamma[i + 12] % 16) == 0);
#endif
		__m128 l1 = _mm_load_ps((float *)&colors->alpha[i + 0]);
		__m128 l2 = _mm_load_ps((float *)&colors->alpha[i + 4]);
		__m128 l3 = _mm_load_ps((float *)&colors->alpha[i + 8]);
		__m128 l4 = _mm_load_ps((float *)&colors->alpha[i + 12]);

		__m128 a1 = _mm_load_ps((float *)&colors->beta[i + 0]);
		__m128 a2 = _mm_load_ps((float *)&colors->beta[i + 4]);
		__m128 a3 = _mm_load_ps((float *)&colors->beta[i + 8]);
		__m128 a4 = _mm_load_ps((float *)&colors->beta[i + 12]);

		__m128 b1 = _mm_load_ps((float *)&colors->gamma[i + 0]);
		__m128 b2 = _mm_load_ps((float *)&colors->gamma[i + 4]);
		__m128 b3 = _mm_load_ps((float *)&colors->gamma[i + 8]);
		__m128 b4 = _mm_load_ps((float *)&colors->gamma[i + 12]);

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
		l += colors->alpha[i];
		a += colors->beta[i];
		b += colors->gamma[i];
	}

	float scale		= 1.0f / num_col;
	cents->alpha[which_cent] = l * scale;
	cents->beta[which_cent]	= a * scale;
	cents->gamma[which_cent] = b * scale;
}

void cielab_avg_sse_cw(const struct Color *restrict colors, struct Color *restrict cents,
		       uint8_t which_cent) {
	uint16_t num_col = colors->num;
	__m128	 sum_l1	 = _mm_setzero_ps();
	__m128	 sum_a1	 = _mm_setzero_ps();
	__m128	 sum_b1	 = _mm_setzero_ps();

	__m128 sum_l2 = _mm_setzero_ps();
	__m128 sum_a2 = _mm_setzero_ps();
	__m128 sum_b2 = _mm_setzero_ps();

	__m128 sum_l3 = _mm_setzero_ps();
	__m128 sum_a3 = _mm_setzero_ps();
	__m128 sum_b3 = _mm_setzero_ps();

	__m128 sum_l4 = _mm_setzero_ps();
	__m128 sum_a4 = _mm_setzero_ps();
	__m128 sum_b4 = _mm_setzero_ps();

	__m128 sum_w = _mm_setzero_ps();

	size_t i = 0;
	for (; i + 15 < num_col; i += 16) {
#ifdef PALETTE_DEBUG
		assert(((uintptr_t)colors->alpha[i + 0] % 32) == 0);
		assert(((uintptr_t)colors->beta[i + 0] % 32) == 0);
		assert(((uintptr_t)colors->gamma[i + 0] % 32) == 0);

		assert(((uintptr_t)colors->alpha[i + 4] % 32) == 0);
		assert(((uintptr_t)colors->beta[i + 4] % 32) == 0);
		assert(((uintptr_t)colors->gamma[i + 4] % 32) == 0);

		assert(((uintptr_t)colors->alpha[i + 8] % 32) == 0);
		assert(((uintptr_t)colors->beta[i + 8] % 32) == 0);
		assert(((uintptr_t)colors->gamma[i + 8] % 32) == 0);

		assert(((uintptr_t)colors->alpha[i + 12] % 32) == 0);
		assert(((uintptr_t)colors->beta[i + 12] % 32) == 0);
		assert(((uintptr_t)colors->gamma[i + 12] % 32) == 0);
#endif
		__m128 l1 = _mm_load_ps(&colors->alpha[i + 0]);
		__m128 l2 = _mm_load_ps(&colors->alpha[i + 4]);
		__m128 l3 = _mm_load_ps(&colors->alpha[i + 8]);
		__m128 l4 = _mm_load_ps(&colors->alpha[i + 12]);

		__m128 a1 = _mm_load_ps(&colors->beta[i + 0]);
		__m128 a2 = _mm_load_ps(&colors->beta[i + 4]);
		__m128 a3 = _mm_load_ps(&colors->beta[i + 8]);
		__m128 a4 = _mm_load_ps(&colors->beta[i + 12]);

		__m128 b1 = _mm_load_ps(&colors->gamma[i + 0]);
		__m128 b2 = _mm_load_ps(&colors->gamma[i + 4]);
		__m128 b3 = _mm_load_ps(&colors->gamma[i + 8]);
		__m128 b4 = _mm_load_ps(&colors->gamma[i + 12]);

		__m128 chroma1 = _mm_sqrt_ps(
		    _mm_add_ps(_mm_mul_ps(a1, a1), _mm_mul_ps(b1, b1)));
		__m128 chroma2 = _mm_sqrt_ps(
		    _mm_add_ps(_mm_mul_ps(a2, a2), _mm_mul_ps(b2, b2)));
		__m128 chroma3 = _mm_sqrt_ps(
		    _mm_add_ps(_mm_mul_ps(a3, a3), _mm_mul_ps(b3, b3)));
		__m128 chroma4 = _mm_sqrt_ps(
		    _mm_add_ps(_mm_mul_ps(a4, a4), _mm_mul_ps(b4, b4)));

		sum_l1 = _mm_add_ps(sum_l1, l1);
		sum_l2 = _mm_add_ps(sum_l2, l2);
		sum_l3 = _mm_add_ps(sum_l3, l3);
		sum_l4 = _mm_add_ps(sum_l4, l4);

		sum_a1 = _mm_add_ps(sum_a1, _mm_mul_ps(a1, chroma1));
		sum_a2 = _mm_add_ps(sum_a2, _mm_mul_ps(a2, chroma2));
		sum_a3 = _mm_add_ps(sum_a3, _mm_mul_ps(a3, chroma3));
		sum_a4 = _mm_add_ps(sum_a4, _mm_mul_ps(a4, chroma4));

		sum_b1 = _mm_add_ps(sum_b1, _mm_mul_ps(b1, chroma1));
		sum_b2 = _mm_add_ps(sum_b2, _mm_mul_ps(b2, chroma2));
		sum_b3 = _mm_add_ps(sum_b3, _mm_mul_ps(b3, chroma3));
		sum_b4 = _mm_add_ps(sum_b4, _mm_mul_ps(b4, chroma4));

		__m128 chroma_sum = _mm_add_ps(_mm_add_ps(chroma1, chroma2),
					       _mm_add_ps(chroma3, chroma4));

		sum_w = _mm_add_ps(sum_w, chroma_sum);
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

	float w_vec[4];
	_mm_store_ps(w_vec, sum_w);
	float w = w_vec[0] + w_vec[1] + w_vec[2] + w_vec[3];

	for (; i < num_col; ++i) {
		float chroma =
		    sqrtf(SQUARE(colors->beta[i]) + SQUARE(colors->gamma[i]));
		l += colors->alpha[i];
		a += colors->beta[i] * chroma;
		b += colors->gamma[i] * chroma;
		w += chroma;
	}

	l /= num_col;
	if (w == 0.0f) {
		a = 0.0f;
		b = 0.0f;
	} else {
		a /= w;
		b /= w;
	}

	cents->alpha[which_cent] = l / num_col;
	cents->beta[which_cent]	= a;
	cents->gamma[which_cent] = b;
}

void oklab_avg_sse(const struct Color *restrict colors, struct Color *restrict cents,
		   uint8_t which_cent) {
	uint16_t num_col = colors->num;
	__m128	 sum_l1	 = _mm_setzero_ps();
	__m128	 sum_a1	 = _mm_setzero_ps();
	__m128	 sum_b1	 = _mm_setzero_ps();

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
		assert(((uintptr_t)colors->alpha[i + 0] % 16) == 0);
		assert(((uintptr_t)colors->beta[i + 0] % 16) == 0);
		assert(((uintptr_t)colors->gamma[i + 0] % 16) == 0);

		assert(((uintptr_t)colors->alpha[i + 4] % 16) == 0);
		assert(((uintptr_t)colors->beta[i + 4] % 16) == 0);
		assert(((uintptr_t)colors->gamma[i + 4] % 16) == 0);

		assert(((uintptr_t)colors->alpha[i + 8] % 16) == 0);
		assert(((uintptr_t)colors->beta[i + 8] % 16) == 0);
		assert(((uintptr_t)colors->gamma[i + 8] % 16) == 0);

		assert(((uintptr_t)colors->alpha[i + 12] % 16) == 0);
		assert(((uintptr_t)colors->beta[i + 12] % 16) == 0);
		assert(((uintptr_t)colors->gamma[i + 12] % 16) == 0);
#endif
		__m128 l1 = _mm_load_ps(&colors->alpha[i + 0]);
		__m128 l2 = _mm_load_ps(&colors->alpha[i + 4]);
		__m128 l3 = _mm_load_ps(&colors->alpha[i + 8]);
		__m128 l4 = _mm_load_ps(&colors->alpha[i + 12]);

		__m128 a1 = _mm_load_ps(&colors->beta[i + 0]);
		__m128 a2 = _mm_load_ps(&colors->beta[i + 4]);
		__m128 a3 = _mm_load_ps(&colors->beta[i + 8]);
		__m128 a4 = _mm_load_ps(&colors->beta[i + 12]);

		__m128 b1 = _mm_load_ps(&colors->gamma[i + 0]);
		__m128 b2 = _mm_load_ps(&colors->gamma[i + 4]);
		__m128 b3 = _mm_load_ps(&colors->gamma[i + 8]);
		__m128 b4 = _mm_load_ps(&colors->gamma[i + 12]);

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
		l += colors->alpha[i];
		a += colors->beta[i];
		b += colors->gamma[i];
	}

	float scale		= 1.0f / num_col;
	cents->alpha[which_cent] = l * scale;
	cents->beta[which_cent]	= a * scale;
	cents->gamma[which_cent] = b * scale;
}

void oklab_avg_sse_cw(const struct Color *restrict colors, struct Color *restrict cents,
		      uint8_t which_cent) {
	uint16_t num_col = colors->num;
	__m128	 sum_l1	 = _mm_setzero_ps();
	__m128	 sum_a1	 = _mm_setzero_ps();
	__m128	 sum_b1	 = _mm_setzero_ps();

	__m128 sum_l2 = _mm_setzero_ps();
	__m128 sum_a2 = _mm_setzero_ps();
	__m128 sum_b2 = _mm_setzero_ps();

	__m128 sum_l3 = _mm_setzero_ps();
	__m128 sum_a3 = _mm_setzero_ps();
	__m128 sum_b3 = _mm_setzero_ps();

	__m128 sum_l4 = _mm_setzero_ps();
	__m128 sum_a4 = _mm_setzero_ps();
	__m128 sum_b4 = _mm_setzero_ps();

	__m128 sum_w = _mm_setzero_ps();

	size_t i = 0;
	for (; i + 15 < num_col; i += 16) {
#ifdef PALETTE_DEBUG
		assert(((uintptr_t)colors->alpha[i + 0] % 32) == 0);
		assert(((uintptr_t)colors->beta[i + 0] % 32) == 0);
		assert(((uintptr_t)colors->gamma[i + 0] % 32) == 0);

		assert(((uintptr_t)colors->alpha[i + 4] % 32) == 0);
		assert(((uintptr_t)colors->beta[i + 4] % 32) == 0);
		assert(((uintptr_t)colors->gamma[i + 4] % 32) == 0);

		assert(((uintptr_t)colors->alpha[i + 8] % 32) == 0);
		assert(((uintptr_t)colors->beta[i + 8] % 32) == 0);
		assert(((uintptr_t)colors->gamma[i + 8] % 32) == 0);

		assert(((uintptr_t)colors->alpha[i + 12] % 32) == 0);
		assert(((uintptr_t)colors->beta[i + 12] % 32) == 0);
		assert(((uintptr_t)colors->gamma[i + 12] % 32) == 0);
#endif
		__m128 l1 = _mm_load_ps(&colors->alpha[i + 0]);
		__m128 l2 = _mm_load_ps(&colors->alpha[i + 4]);
		__m128 l3 = _mm_load_ps(&colors->alpha[i + 8]);
		__m128 l4 = _mm_load_ps(&colors->alpha[i + 12]);

		__m128 a1 = _mm_load_ps(&colors->beta[i + 0]);
		__m128 a2 = _mm_load_ps(&colors->beta[i + 4]);
		__m128 a3 = _mm_load_ps(&colors->beta[i + 8]);
		__m128 a4 = _mm_load_ps(&colors->beta[i + 12]);

		__m128 b1 = _mm_load_ps(&colors->gamma[i + 0]);
		__m128 b2 = _mm_load_ps(&colors->gamma[i + 4]);
		__m128 b3 = _mm_load_ps(&colors->gamma[i + 8]);
		__m128 b4 = _mm_load_ps(&colors->gamma[i + 12]);

		__m128 chroma1 = _mm_sqrt_ps(
		    _mm_add_ps(_mm_mul_ps(a1, a1), _mm_mul_ps(b1, b1)));
		__m128 chroma2 = _mm_sqrt_ps(
		    _mm_add_ps(_mm_mul_ps(a2, a2), _mm_mul_ps(b2, b2)));
		__m128 chroma3 = _mm_sqrt_ps(
		    _mm_add_ps(_mm_mul_ps(a3, a3), _mm_mul_ps(b3, b3)));
		__m128 chroma4 = _mm_sqrt_ps(
		    _mm_add_ps(_mm_mul_ps(a4, a4), _mm_mul_ps(b4, b4)));

		sum_l1 = _mm_add_ps(sum_l1, l1);
		sum_l2 = _mm_add_ps(sum_l2, l2);
		sum_l3 = _mm_add_ps(sum_l3, l3);
		sum_l4 = _mm_add_ps(sum_l4, l4);

		sum_a1 = _mm_add_ps(sum_a1, _mm_mul_ps(a1, chroma1));
		sum_a2 = _mm_add_ps(sum_a2, _mm_mul_ps(a2, chroma2));
		sum_a3 = _mm_add_ps(sum_a3, _mm_mul_ps(a3, chroma3));
		sum_a4 = _mm_add_ps(sum_a4, _mm_mul_ps(a4, chroma4));

		sum_b1 = _mm_add_ps(sum_b1, _mm_mul_ps(b1, chroma1));
		sum_b2 = _mm_add_ps(sum_b2, _mm_mul_ps(b2, chroma2));
		sum_b3 = _mm_add_ps(sum_b3, _mm_mul_ps(b3, chroma3));
		sum_b4 = _mm_add_ps(sum_b4, _mm_mul_ps(b4, chroma4));

		__m128 chroma_sum = _mm_add_ps(_mm_add_ps(chroma1, chroma2),
					       _mm_add_ps(chroma3, chroma4));

		sum_w = _mm_add_ps(sum_w, chroma_sum);
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

	float w_vec[4];
	_mm_store_ps(w_vec, sum_w);
	float w = w_vec[0] + w_vec[1] + w_vec[2] + w_vec[3];

	for (; i < num_col; ++i) {
		float chroma =
		    sqrtf(SQUARE(colors->beta[i]) + SQUARE(colors->gamma[i]));
		l += colors->alpha[i];
		a += colors->beta[i] * chroma;
		b += colors->gamma[i] * chroma;
		w += chroma;
	}

	l /= num_col;
	if (w == 0.0f) {
		a = 0.0f;
		b = 0.0f;
	} else {
		a /= w;
		b /= w;
	}

	cents->alpha[which_cent] = l / num_col;
	cents->beta[which_cent]	= a;
	cents->gamma[which_cent] = b;
}
