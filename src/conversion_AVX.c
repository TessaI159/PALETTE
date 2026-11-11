#include "Color.h"
#include "Conversion_Internal.h"

#if (defined(__x86_64__) || defined(_M_X64) || defined(__i386)) &&             \
    defined(__AVX)
    #include <immintrin.h>

/* Default to SSE even if AVX is available for now to resolve linker
 * issues */
void _mm256_srgb_to_oklab_ps(struct Color *restrict colors) {
  cvt_mm_srgb_to_oklab_ps(colors);
}

void _mm256_oklab_to_srgb_ps(struct Color *restrict colors) {
  cvt_mm_oklab_to_srgb_ps(colors);
}

/* Default to fallbacks even if AVX is available for now to resolve linker
 * issues */
void _mm256_srgb_to_cielab_ps(struct Color *restrict colors) {
  cielab_to_srgb_fb(colors);
}
void _mm256_cielab_to_srgb_ps(struct Color *restrict colors) {
  cielab_to_srgb_fb(colors);
}

#else

void _mm256_srgb_to_oklab_ps(struct Color *restrict colors) {
  srgb_to_oklab_fb(colors);
}

void _mm256_oklab_to_srgb_ps(struct Color *restrict colors) {
  oklab_to_srgb_fb(colors);
}

void _mm256_srgb_to_cielab_ps(struct Color *restrict colors) {
  cielab_to_srgb_fb(colors);
}
void _mm256_cielab_to_srgb_ps(struct Color *restrict colors) {
  cielab_to_srgb_fb(colors);
}

#endif
