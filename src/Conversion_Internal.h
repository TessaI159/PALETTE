#ifndef CONVERSION_INTERNAL_H
#define CONVERSION_INTERNAL_H

typedef struct Color Color;

void srgb_to_oklab_fb(struct Color *restrict colors);
void srgb_to_cielab_fb(struct Color *restrict colors);
void oklab_to_srgb_fb(struct Color *restrict colors);
void cielab_to_srgb_fb(struct Color *restrict colors);

void _mm_srgb_to_oklab_ps(struct Color *restrict colors);
void _mm_srgb_to_cielab_ps(struct Color *restrict colors);
void _mm_oklab_to_srgb_ps(struct Color *restrict colors);
void _mm_cielab_to_srgb_ps(struct Color *restrict colors);

void _mm256_srgb_to_oklab_ps(struct Color *restrict colors);
void _mm256_srgb_to_cielab_ps(struct Color *restrict colors);
void _mm256_oklab_to_srgb_ps(struct Color *restrict colors);
void _mm256_cielab_to_srgb_ps(struct Color *restrict colors);

#endif
