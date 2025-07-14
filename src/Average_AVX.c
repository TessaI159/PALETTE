#ifdef PALETTE_DEBUG
#include <assert.h>
#endif
#include <stdio.h>
#ifdef __AVX__
#include <immintrin.h>
#endif

#include "Color.h"

struct Color cielab_avg_avx(const struct cielab_SoA *colors, uint16_t num_col) {

	__m256	 lsum = _mm256_setzero_ps();
	__m256	 asum = _mm256_setzero_ps();
	__m256	 bsum = _mm256_setzero_ps();
	uint16_t i;
	for (i = 0; i + 7 < num_col; i += 8) {
#ifdef PALETTE_DEBUG
		assert(((uintptr_t)&colors->l[i] % 32) == 0);
		assert(((uintptr_t)&colors->a[i] % 32) == 0);
		assert(((uintptr_t)&colors->b[i] % 32) == 0);
#endif
		__m256 l = _mm256_load_ps(&colors->l[i]);
		__m256 a = _mm256_load_ps(&colors->a[i]);
		__m256 b = _mm256_load_ps(&colors->b[i]);
	}
	return Color_create(0, 0.0f, 0.0f);
}

/* struct Color cielab_avg_avx(struct Color *colors, uint16_t num_col) { */
/* 	__m256 lsum = _mm256_setzero_ps(); */
/* 	__m256 asum = _mm256_setzero_ps(); */
/* 	__m256 bsum = _mm256_setzero_ps(); */

/* 	uint16_t i; */
/* 	for (i = 0; i + 8 <= num_col; i += 8) { */
/* 		__m256 l = _mm256_set_ps( */
/* 		    colors[i + 7].data.cielab.l, colors[i + 6].data.cielab.l, */
/* 		    colors[i + 5].data.cielab.l, colors[i + 4].data.cielab.l, */
/* 		    colors[i + 3].data.cielab.l, colors[i + 2].data.cielab.l, */
/* 		    colors[i + 1].data.cielab.l, colors[i + 0].data.cielab.l);
 */
/* 		__m256 a = _mm256_set_ps( */
/* 		    colors[i + 7].data.cielab.a, colors[i + 6].data.cielab.a, */
/* 		    colors[i + 5].data.cielab.a, colors[i + 4].data.cielab.a, */
/* 		    colors[i + 3].data.cielab.a, colors[i + 2].data.cielab.a, */
/* 		    colors[i + 1].data.cielab.a, colors[i + 0].data.cielab.a);
 */
/* 		__m256 b = _mm256_set_ps( */
/* 		    colors[i + 7].data.cielab.b, colors[i + 6].data.cielab.b, */
/* 		    colors[i + 5].data.cielab.b, colors[i + 4].data.cielab.b, */
/* 		    colors[i + 3].data.cielab.b, colors[i + 2].data.cielab.b, */
/* 		    colors[i + 1].data.cielab.b, colors[i + 0].data.cielab.b);
 */

/* 		lsum = _mm256_add_ps(lsum, l); */
/* 		asum = _mm256_add_ps(asum, a); */
/* 		bsum = _mm256_add_ps(bsum, b); */
/* 	} */

/* 	// Horizontal sum of AVX registers */
/* 	__m128 llow    = _mm256_castps256_ps128(lsum); */
/* 	__m128 lhigh   = _mm256_extractf128_ps(lsum, 1); */
/* 	__m128 lsum128 = _mm_add_ps(llow, lhigh); */
/* 	lsum128	       = _mm_hadd_ps(lsum128, lsum128); */
/* 	lsum128	       = _mm_hadd_ps(lsum128, lsum128); */
/* 	float l_total  = _mm_cvtss_f32(lsum128); */

/* 	__m128 alow    = _mm256_castps256_ps128(asum); */
/* 	__m128 ahigh   = _mm256_extractf128_ps(asum, 1); */
/* 	__m128 asum128 = _mm_add_ps(alow, ahigh); */
/* 	asum128	       = _mm_hadd_ps(asum128, asum128); */
/* 	asum128	       = _mm_hadd_ps(asum128, asum128); */
/* 	float a_total  = _mm_cvtss_f32(asum128); */

/* 	__m128 blow    = _mm256_castps256_ps128(bsum); */
/* 	__m128 bhigh   = _mm256_extractf128_ps(bsum, 1); */
/* 	__m128 bsum128 = _mm_add_ps(blow, bhigh); */
/* 	bsum128	       = _mm_hadd_ps(bsum128, bsum128); */
/* 	bsum128	       = _mm_hadd_ps(bsum128, bsum128); */
/* 	float b_total  = _mm_cvtss_f32(bsum128); */

/* 	for (; i < num_col; ++i) { */
/* 		l_total += colors[i].data.cielab.l; */
/* 		a_total += colors[i].data.cielab.a; */
/* 		b_total += colors[i].data.cielab.b; */
/* 	} */

/* 	return Color_create_cielab(l_total / (float)num_col, */
/* 				   a_total / (float)num_col, */
/* 				   b_total / (float)num_col); */
/* } */
