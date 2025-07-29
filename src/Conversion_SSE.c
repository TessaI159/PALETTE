#include <math.h>

#include "Color.h"
#include "Math_SSE.h"

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

static inline void srgb_to_oklab_scalar(struct Color *restrict colors,
					int start_index) {
	const uint32_t num_col = colors->num;
	for (uint32_t i = start_index; i < num_col; ++i) {
		double r = colors->alpha[i];
		double g = colors->beta[i];
		double b = colors->gamma[i];
		double l =
		    0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b;
		double m =
		    0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b;
		double s =
		    0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b;
		l = cbrt(l);
		m = cbrt(m);
		s = cbrt(s);
		/* cbrt has domain [0, 1] */

		colors->alpha[i] =
		    0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s;
		colors->beta[i] =
		    1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s;
		colors->gamma[i] =
		    0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s;
	}
}

void _mm_srgb_to_oklab_ps(struct Color *restrict colors) {
	const uint32_t num_col = colors->num;
	const __m128   zero    = _mm_set1_ps(0.0f);
	uint32_t       i       = 0;
	for (; i + FL_PER_VEC <= num_col; i += FL_PER_VEC) {
		__m128 r = _mm_loadu_ps(&colors->alpha[i]);
		__m128 g = _mm_loadu_ps(&colors->beta[i]);
		__m128 b = _mm_loadu_ps(&colors->gamma[i]);
		__m128 l =
		    FMADD(_mm_set1_ps(0.4122214708f), r,
			  FMADD(_mm_set1_ps(0.5363325363f), g,
				FMADD(_mm_set1_ps(0.0514459929f), b, zero)));
		__m128 m =
		    FMADD(_mm_set1_ps(0.2119034982f), r,
			  FMADD(_mm_set1_ps(0.6806995451f), g,
				FMADD(_mm_set1_ps(0.1073969566f), b, zero)));
		__m128 s =
		    FMADD(_mm_set1_ps(0.0883024619f), r,
			  FMADD(_mm_set1_ps(0.2817188376f), g,
				FMADD(_mm_set1_ps(0.6299787005f), b, zero)));

		l = _mm_cbrt_ps(l);
		m = _mm_cbrt_ps(m);
		s = _mm_cbrt_ps(s);

		r = FMADD(_mm_set1_ps(0.2104542553f), l,
			  FMADD(_mm_set1_ps(0.7936177850f), m,
				FMADD(_mm_set1_ps(-0.0040720468f), s, zero)));
		g = FMADD(_mm_set1_ps(1.9779984951f), l,
			  FMADD(_mm_set1_ps(-2.4285922050f), m,
				FMADD(_mm_set1_ps(0.4505937099f), s, zero)));
		b = FMADD(_mm_set1_ps(0.0259040371f), l,
			  FMADD(_mm_set1_ps(0.7827717662f), m,
				FMADD(_mm_set1_ps(-0.8086757660f), s, zero)));

		_mm_storeu_ps(&colors->alpha[i], r);
		_mm_storeu_ps(&colors->beta[i], g);
		_mm_storeu_ps(&colors->gamma[i], b);
	}
	srgb_to_oklab_scalar(colors, i);
}
