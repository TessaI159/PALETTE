#ifndef CONVERSION_H
#define CONVERSION_H

struct Color Color;
/* Conversions between linear normalized sRGB [0,1] and perceptual spaces */
void srgb_to_oklab(struct Color *restrict colors);
void srgb_to_cielab(struct Color *restrict colors);
void oklab_to_srgb(struct Color *restrict colors);
void cielab_to_srgb(struct Color *restrict colors);

//* Gamma encode and decode. [0, 1] */
void srgb_linearize(struct Color *restrict colors);
void srgb_gamma_correct(struct Color *restrict colors);

#endif
