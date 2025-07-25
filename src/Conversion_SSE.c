#include "Math_SSE.h"

#include "Color.h"

#define V(x) _mm_set1_ps((x))
#define FMADD(a, b, c) _mm_add_ps(_mm_mul_ps((a), (b)), (c))
#define FL_PER_VEC 4

static inline void _mm_srgb_to_xyz_ps(const __m128 *restrict r,
				      const __m128 *restrict g,
				      const __m128 *restrict b,
				      __m128 *restrict x, __m128 *restrict y,
				      __m128 *restrict z) {
	*x = FMADD(*r, V(0.4124f),
		   FMADD(*g, V(0.3576f), _mm_mul_ps(*b, V(0.1805f))));
	*y = FMADD(*r, V(0.2126f),
		   FMADD(*g, V(0.7152f), _mm_mul_ps(*b, V(0.0722f))));
	*z = FMADD(*r, V(0.0193f),
		   FMADD(*g, V(0.1192f), _mm_mul_ps(*b, V(0.9505f))));
}

static inline void _mm_xyz_to_srgb_ps(const __m128 *restrict x,
				      const __m128 *restrict y,
				      const __m128 *restrict z,
				      __m128 *restrict r, __m128 *restrict g,
				      __m128 *restrict b) {
	*r = FMADD(*x, V(3.2406f),
		   FMADD(*y, V(-1.5372f), _mm_mul_ps(*z, V(-0.4986f))));
	*g = FMADD(*x, V(-0.9689f),
		   FMADD(*y, V(1.8758f), _mm_mul_ps(*z, V(0.0415f))));
	*b = FMADD(*x, V(0.0557f),
		   FMADD(*y, V(-0.2040f), _mm_mul_ps(*z, V(1.0570f))));
}

void _mm_srgb_to_oklab_ps(struct Color *restrict colors) {
	const uint32_t num_col = colors->num;
	for (uint32_t i = 0; i + FL_PER_VEC <= num_col; i += FL_PER_VEC) {
		__m128 r = _mm_loadu_ps(colors->alpha);
		__m128 g = _mm_loadu_ps(colors->beta);
		__m128 b = _mm_loadu_ps(colors->gamma);
	}
}
