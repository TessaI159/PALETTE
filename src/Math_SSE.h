#ifndef MATH_SSE_H
#define MATH_SSE_H
#include <pmmintrin.h>
#include <smmintrin.h>

#define _mm_abs_ps(x) _mm_and_ps((x), _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF)))
#define FMADD(a, b, c) _mm_add_ps(_mm_mul_ps((a), (b)), (c))
#define _MM_PI 3.14159265f
#define _MM_2PI 6.28318531f
#define _MM_PI_2 1.57079633f
#define _MM_1_OVER_2PI 0.159154943f

__m128 _mm_pow7_ps(__m128 x);
__m128 _mm_square_ps(__m128 x);
__m128 _mm_deg_2_rad_ps(__m128 x);
__m128 _mm_rad_2_deg_ps(__m128 x);
__m128 _mm_atan_ps(__m128 x);
__m128 _mm_atan2_ps(__m128 y, __m128 x);
__m128 _mm_sin_ps(__m128 x);
__m128 _mm_cos_ps(__m128 x);
__m128 _mm_exp_ps(__m128 x);
__m128 _mm_cbrt_ps(__m128 x);

#endif
