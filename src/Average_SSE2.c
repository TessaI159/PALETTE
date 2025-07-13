#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#if defined(__x86_64__) || defined(_M_X64)
#include <emmintrin.h>
#endif

#include "Color.h"

extern void convert_to(struct Color *color, enum ColorSpace);
typedef __m128 (*_load)(const float *);

static inline __m128 load_aligned(const float *p) {
	return _mm_load_ps(p);
}

static inline __m128 load_unaligned(const float *p) {
	return _mm_loadu_ps(p);
}

struct Color cielab_avg_sse2(struct Color *colors, uint16_t num_col) {
	__m128 sum = _mm_setzero_ps();

	_load load;
	if (offsetof(struct Color, data) == 0) {
	  load = load_aligned;
	} else {
	  load = load_unaligned;
	}
	
	size_t i = 0;
	for (; i + 3 < num_col; i += 4) {
		__m128 lab1 =
		    load((const float *)&colors[i + 0].data.cielab);
		__m128 lab2 =
		    load((const float *)&colors[i + 1].data.cielab);
		__m128 lab3 =
		    load((const float *)&colors[i + 2].data.cielab);
		__m128 lab4 =
		    load((const float *)&colors[i + 3].data.cielab);

		__m128 tmp1 = _mm_add_ps(lab1, lab2);
		__m128 tmp2 = _mm_add_ps(lab3, lab4);
		sum	    = _mm_add_ps(sum, _mm_add_ps(tmp1, tmp2));
	}
	for (; i < num_col; ++i) {
		__m128 lab = load((float *)&colors[i].data.cielab);
		sum	   = _mm_add_ps(sum, lab);
	}
	float accum[4];
	_mm_store_ps(accum, sum);

	float l = accum[0] / num_col;
	float a = accum[1] / num_col;
	float b = accum[2] / num_col;

	return Color_create(l, a, b);
	return Color_create(0, 0, 0);
}
