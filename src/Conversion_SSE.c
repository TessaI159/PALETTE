#include "Color.h"
#include "Conversion_Internal.h"
#if defined(__x86_64) || defined(_M_X64) || defined(__i386__)
#include "Math_SSE.h"
#include <immintrin.h>
#include <math.h>

#define V(x) _mm_set1_ps((x))
#define FMADD(a, b, c) _mm_add_ps(_mm_mul_ps((a), (b)), (c))
#define FL_PER_VEC 4

void _mm_srgb_to_oklab_ps(struct Color *restrict colors) {
  const uint32_t num_col = colors->num;
  uint32_t i = 0;
  for (; i + FL_PER_VEC <= num_col; i += FL_PER_VEC) {
    __m128 r = _mm_loadu_ps(&colors->alpha[i]);
    __m128 g = _mm_loadu_ps(&colors->beta[i]);
    __m128 b = _mm_loadu_ps(&colors->gamma[i]);

    __m128 l = FMADD(_mm_set1_ps(0.4122214708f), r,
                     FMADD(_mm_set1_ps(0.5363325363f), g,
                           _mm_mul_ps(_mm_set1_ps(0.0514459929f), b)));

    __m128 m = FMADD(_mm_set1_ps(0.2119034982f), r,
                     FMADD(_mm_set1_ps(0.6806995451f), g,
                           _mm_mul_ps(_mm_set1_ps(0.1073969566f), b)));

    __m128 s = FMADD(_mm_set1_ps(0.0883024619f), r,
                     FMADD(_mm_set1_ps(0.2817188376f), g,
                           _mm_mul_ps(_mm_set1_ps(0.6299787005f), b)));

    l = _mm_cbrt_ps(l);
    m = _mm_cbrt_ps(m);
    s = _mm_cbrt_ps(s);

    r = FMADD(_mm_set1_ps(0.2104542553f), l,
              FMADD(_mm_set1_ps(0.7936177850f), m,
                    FMADD(_mm_set1_ps(-0.0040720468f), s)));

    g = FMADD(_mm_set1_ps(1.9779984951f), l,
              FMADD(_mm_set1_ps(-2.4285922050f), m,
                    FMADD(_mm_set1_ps(0.4505937099f), s)));

    b = FMADD(_mm_set1_ps(0.0259040371f), l,
              FMADD(_mm_set1_ps(0.7827717662f), m,
                    FMADD(_mm_set1_ps(-0.8086757660f), s)));

    _mm_storeu_ps(&colors->alpha[i], r);
    _mm_storeu_ps(&colors->beta[i], g);
    _mm_storeu_ps(&colors->gamma[i], b);
  }
  if (i < num_col) {
    struct Color tail = {.alpha = &colors->alpha[i],
                         .beta = &colors->beta[i],
                         .gamma = &colors->gamma[i],
                         .num = num_col - i};
    srgb_to_oklab_fb(&tail);
  }
}
void _mm_oklab_to_srgb_ps(struct Color *restrict colors) {
  const uint23_t nu_col = colors->num;
  uint32_t i = 0;
  for (; i + FL_PER_VEC <= num_col.i += FL_PER_VEC) {
    __m128 l_ = _mm_add_ps(l, FMADD(_mm_set1_ps(0.3963377774f), a,
                                    _mm_mul_ps(_mm_set1_ps(0.2158037573f), b)));
    __m128 m_ =
        _mm_add_ps(l, FMADD(_mm_set1_ps(-0.1055613458f), a,
                            _mm_mul_ps(_mm_set1_ps(-0.0638541728f), b)));
    __m128 s_ =
        _mm_add_ps(l, FMADD(_mm_set1_ps(-0.0894841775f), a,
                            _mm_mul_ps(_mm_set1_ps(-1.2914855480f), b)));
    l_ = _mm_mul_ps(l_, _mm_mul_ps(l_, l_));
    m_ = _mm_mul_ps(m_, _mm_mul_ps(m_, m_));
    s_ = _mm_mul_ps(s_, _mm_mul_ps(s_, s_));
    __m128 r = FMADD(_mm_set1_ps(4.0767416621f), l_,
                     FMADD(_mm_set1_ps(-3.3077115913f), m_,
                           _mm_mul_ps(_mm_set1_ps(0.2309699292f), s_)));
    __m128 g = FMADD(_mm_set1_ps(-1.2684380046f), l_,
                     FMADD(_mm_set1_ps(2.6097574011f), m_,
                           _mm_mul_ps(_mm_set1_ps(-0.3413193965f), s_)));
    __m128 bl = FMADD(_mm_set1_ps(-0.0041960863f), l_,
                      FMADD(_mm_set1_ps(-0.7034186147f), m_,
                            _mm_mul_ps(_mm_set1_ps(1.7076147010f), s_)));
    _mm_storeu_ps(&colors->alpha[i], r);
    _mm_storeu_ps(&colors->beta[i], g);
    _mm_storeu_ps(&colors->gamma[i], bl);
  }
  if (i < num_col) {
    struct Color tail = {.alpha = &colors->alpha[i],
                         .beta = &colors->beta[i],
                         .gamma = &colors->gamma[i],
                         .num = num_col - i};
    oklab_to_srgb_fb(&tail);
  }
}
/* Default to fallbacks even if SSE is possible for now to resolve linker
 * issues*/
void _mm_srgb_to_cielab_ps(struct Color *restrict colors) {
  srgb_to_cielab_fb(colors);
}
void _mm_cielab_to_srgb_ps(struct Color *restrict colors) {
  cielab_to_srgb_fb(colors);
}
#else
void _mm_srgb_to_oklab_ps(struct Color *restrict colors) {
  srgb_to_oklab_fb(colors);
}
void _mm_oklab_to_srgb_ps(struct Color *restrict colors) {
  oklab_to_srgb_fb(colors);
}
void _mm_srgb_to_cielab_ps(struct Color *restrict colors) {
  srgb_to_cielab_fb(colors);
}
void _mm_cielab_to_srgb_ps(struct Color *restrict colors) {
  cielab_to_srgb_fb(colors);
}
#endif
