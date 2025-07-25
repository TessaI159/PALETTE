#ifndef CONVERSION_H
#define CONVERSION_H

typedef struct Color Color;

void srgb_to_oklab(struct Color *restrict colors);
void srgb_to_cielab(struct Color *restrict colors);
void oklab_to_srgb(struct Color *restrict colors);
void cielab_to_srgb(struct Color *restrict colors);

#endif
